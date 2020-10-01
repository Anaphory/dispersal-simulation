import collections
import typing as t
from pathlib import Path

import json
import more_itertools
import numpy

import rasterio
import rasterio.features
from h3.api import numpy_int as h3
import zipfile

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep
from shapely.validation import make_valid

from ecoregions import ECOREGIONS, TC, RESOLUTION

import cartopy.geodesic as geodesic

GEODESIC: geodesic.Geodesic = geodesic.Geodesic()


from raster_data import *
from raster_data import (
    ecoregion_tile_from_geocoordinates,
    gmted_tile_from_geocoordinates,
)


def navigation_speed(slope: float) -> float:
    """Using the slope in %, calculate the navigation speed in m/s

    This function calculates the off-road navigation speed (for male cadets in
    forested areas with navigational targets every 20 minutes) following
    [@irmischer2018measuring]. Like their formula, slope is in % and speed in
    m/s.

    > [T]he fastest off-road navigation speed was 0.78 m/s for males […] with a
    > peak at −2%.

    >>> navigation_speed(-2.)
    0.78
    >>> navigation_speed(-2.) > navigation_speed(-1.5)
    True
    >>> navigation_speed(-2.) > navigation_speed(-2.5)
    True

    """
    return 0.11 + 0.67 * numpy.exp(-((slope + 2.0) ** 2) / 1800.0)


def all_pairwise_distances(
    elevation: numpy.array,
    transform: rasterio.Affine,
    terrain_coefficients: numpy.array,
):
    d_n, d_e, d_ne = [], [], []
    for y in range(len(elevation)):
        (lon0, lat0) = transform * (0, y + 1)
        (lon1, lat1) = transform * (1, y)

        d = GEODESIC.inverse((lon0, lat0), [(lon0, lat1), (lon1, lat0), (lon1, lat1)])
        d_n.append(d[0, 0])
        d_e.append(d[1, 0])
        d_ne.append(d[2, 0])
    distance_to_north = numpy.array(d_n)[:-1]
    slope_to_north = (
        100 * (elevation[1:, :] - elevation[:-1, :]) / distance_to_north[:, None]
    )
    tc_to_north = (terrain_coefficients[1:, :] + terrain_coefficients[:-1, :]) / 2
    north = numpy.full(elevation.shape, numpy.nan)
    north[1:, :] = distance_to_north[:, None]
    north[1:, :] /= navigation_speed(slope_to_north) * tc_to_north
    south = numpy.full(elevation.shape, numpy.nan)
    south[:-1, :] = distance_to_north[:, None]
    south[:-1, :] /= navigation_speed(-slope_to_north) * tc_to_north
    del distance_to_north, slope_to_north, tc_to_north

    distance_to_east = numpy.array(d_e)
    slope_to_east = (
        100 * (elevation[:, 1:] - elevation[:, :-1]) / distance_to_east[:, None]
    )
    tc_to_east = (terrain_coefficients[:, 1:] + terrain_coefficients[:, :-1]) / 2
    east = numpy.full(elevation.shape, numpy.nan)
    east[:, :-1] = distance_to_east[:, None]
    east[:, :-1] /= navigation_speed(slope_to_east) * tc_to_east
    west = numpy.full(elevation.shape, numpy.nan)
    west[:, 1:] = distance_to_east[:, None]
    west[:, 1:] /= navigation_speed(-slope_to_east) * tc_to_east
    del distance_to_east, slope_to_east, tc_to_east

    areas = [n * e for n, e in zip(d_n, d_e)]

    distance_to_northeast = numpy.array(d_ne)[:-1]
    slope_to_northeast = (
        100 * (elevation[1:, 1:] - elevation[:-1, :-1]) / distance_to_northeast[:, None]
    )
    tc_to_northeast = (
        terrain_coefficients[1:, 1:] + terrain_coefficients[:-1, :-1]
    ) / 2
    northeast = numpy.full(elevation.shape, numpy.nan)
    northeast[1:, :-1] = distance_to_northeast[:, None]
    northeast[1:, :-1] /= navigation_speed(slope_to_northeast) * tc_to_northeast
    southwest = numpy.full(elevation.shape, numpy.nan)
    southwest[:-1, 1:] = distance_to_northeast[:, None]
    southwest[:-1, 1:] /= navigation_speed(-slope_to_northeast) * tc_to_northeast
    del distance_to_northeast, slope_to_northeast, tc_to_northeast
    distance_to_northwest = numpy.array(d_ne)[:-1]
    slope_to_northwest = (
        100 * (elevation[1:, :-1] - elevation[:-1, 1:]) / distance_to_northwest[:, None]
    )

    tc_to_northwest = (
        terrain_coefficients[1:, :-1] + terrain_coefficients[:-1, 1:]
    ) / 2

    southeast = numpy.full(elevation.shape, numpy.nan)
    southeast[:-1, :-1] = distance_to_northwest[:, None]
    southeast[:-1, :-1] /= navigation_speed(-slope_to_northwest) * tc_to_northwest
    northwest = numpy.full(elevation.shape, numpy.nan)
    northwest[1:, 1:] = distance_to_northwest[:, None]
    northwest[1:, 1:] /= navigation_speed(slope_to_northwest) * tc_to_northwest
    del distance_to_northwest, slope_to_northwest, tc_to_northwest

    return {
        (-1, 0): north,
        (-1, 1): northeast,
        (0, 1): east,
        (1, 1): southeast,
        (1, 0): south,
        (1, -1): southwest,
        (0, -1): west,
        (-1, -1): northwest,
    }, areas


def travel_time_raster(lon: float, lat: float, bf: int = 500):  # Buffer
    """Compute pairwise distances and ecoregion composition

    Given a digital elevation model as raster map and a matching ecoregions
    raster map, compute the pairwise distance to its 1-hex and 2-hex neighbors
    for every H3 address hex at standard resolution, as well as the approximate
    cover of that cell in terms of ecoregions, for all cells where that is
    possible.

    Distances are computed in gross hours of travel while navigating off-track,
    following [@irmischer2018measuring].

    Returns
    =======
    d: A mapping. d[h1][h2] is the distance, in hours, from the center of h1 to
        the center of h2.
    e: A mapping. d[h1][b] is the proportion of hex h1 covered by ecoregion b.

    """
    elevation_file = gmted_tile_from_geocoordinates(lon, lat)
    m_trafo = elevation_file.transform
    height, width = elevation_file.shape
    elevation = numpy.full((bf + height + bf, bf + width + bf), -100, int)
    ecoregions = numpy.full((bf + height + bf, bf + width + bf), 999, int)
    for delta_lon, target_cols, source_cols in [
        (-30, slice(0, bf), slice(width - bf, width)),
        (0, slice(bf, width + bf), slice(0, width)),
        (+30, slice(width + bf, width + 2 * bf), slice(0, bf)),
    ]:
        for delta_lat, target_rows, source_rows in [
            (+20, slice(0, bf), slice(height - bf, height)),
            (0, slice(bf, height + bf), slice(0, height)),
            (-20, slice(height + bf, height + 2 * bf), slice(0, bf)),
        ]:
            try:
                elevation[target_rows, target_cols] = (
                    gmted_tile_from_geocoordinates(lon + delta_lon, lat + delta_lat)
                ).read(1)[source_rows, source_cols]
            except rasterio.RasterioIOError:
                pass
            try:
                ecoregions[target_rows, target_cols] = (
                    ecoregion_tile_from_geocoordinates(lon + delta_lon, lat + delta_lat)
                ).read(1)[source_rows, source_cols]
            except rasterio.RasterioIOError:
                pass

    print("Computing hex extents…")
    transform = rasterio.Affine(
        m_trafo.a,
        0,
        m_trafo.c - bf * m_trafo.a,
        0,
        m_trafo.e,
        m_trafo.f - bf * m_trafo.e,
    )

    print("Computing terrain coefficients…")
    terrain_coefficient_raster = TC[ecoregions]

    print("Computing distances on the grid…")
    distance_by_direction = all_pairwise_distances(
        elevation, transform, terrain_coefficient_raster
    )

    return distance_by_direction


class RiverNetwork:
    cache = None

    @classmethod
    def reaches(cls):
        """
        >>> rivers = RiverNetwork.reaches()
        >>> rivers.fields
        [('DeletionFlag', 'C', 1, 0), ['Reach_ID', 'N', 10, 0], ['Next_down', 'N', 10, 0], ['Length_km', 'F', 13, 11], ['Log_Q_avg', 'F', 13, 11], ['Log_Q_var', 'F', 13, 11], ['Class_hydr', 'N', 10, 0], ['Temp_min', 'F', 13, 11], ['CMI_indx', 'F', 13, 11], ['Log_elev', 'F', 13, 11], ['Class_phys', 'N', 10, 0], ['Lake_wet', 'N', 10, 0], ['Stream_pow', 'F', 13, 11], ['Class_geom', 'N', 10, 0], ['Reach_type', 'N', 10, 0], ['Kmeans_30', 'N', 10, 0]]
        >>> rivers.numRecords
        ...

        """
        if cls.cache is not None:
            return cls.cache
        zipshape = zipfile.ZipFile(
            (Path(__file__).parent / "../rivers/GloRiC_v10_shapefile.zip").open("rb")
        )
        shape = shapefile.Reader(
            shp=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.shp"),
            shx=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.shx"),
            dbf=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.dbf"),
            encoding="utf-8",
        )
        cls.cache = shape
        return shape

    def __init__(self, mask: t.Optional[sgeom.Polygon] = None):
        """

        >>> eco = RiverNetwork()
        >>> eco.record(508)[1]
        'Tocantins/Pindare moist forests'

        """
        self.shp = self.reaches()


RIVERS = RiverNetwork()

oceans_zip = zipfile.ZipFile(
    (Path(__file__).parent / "../naturalearth/ne_10m_ocean.zip").open("rb")
)
OCEANS = shapefile.Reader(
    shp=oceans_zip.open("ne_10m_ocean.shp"),
    shx=oceans_zip.open("ne_10m_ocean.shx"),
    dbf=oceans_zip.open("ne_10m_ocean.dbf"),
    encoding="utf-8",
)
OCEANS = next(OCEANS.iterShapes())
OCEANS = make_valid(sgeom.shape(OCEANS))

lakes_zip = zipfile.ZipFile(
    (Path(__file__).parent / "../naturalearth/ne_10m_lakes.zip").open("rb")
)
LAKES = shapefile.Reader(
    shp=lakes_zip.open("ne_10m_lakes.shp"),
    shx=lakes_zip.open("ne_10m_lakes.shx"),
    dbf=lakes_zip.open("ne_10m_lakes.dbf"),
    encoding="utf-8",
)

rivers_zip = zipfile.ZipFile(
    (Path(__file__).parent / "../naturalearth/ne_10m_rivers_lake_centerlines.zip").open(
        "rb"
    )
)
MORE_RIVERS = shapefile.Reader(
    shp=rivers_zip.open("ne_10m_rivers_lake_centerlines.shp"),
    shx=rivers_zip.open("ne_10m_rivers_lake_centerlines.shx"),
    dbf=rivers_zip.open("ne_10m_rivers_lake_centerlines.dbf"),
    encoding="utf-8",
)


# Somewhere, I found speeds of 4.5 knots for kayak cruising. That's 8.334 km/h,
# but the database stores data in seconds, so 2.315 m/s. That's about a factor
# 3 faster than walking slightly downhill (before terrain coefficients), which
# sounds sensible.
KAYAK_SPEED = 2.315


def estimate_flow_speed(discharge, slope):
    """Estimate the flow speed, in km/s from discharge and slope

    This is a very rough estimate, following [@schulze2005simulating]. They
    suggest to at least estimate the widely varying river roughness n, but we
    cannot do that, so we use their reported mean value of 0.044.

    W = 2.71 · Q^0.557
    D = 0.349 · Q^0.341
    R = D · W / (2 D + W)
    n = 0.044
    v = 1/n · R^2/3 · S^1/2

    """
    n = 0.044
    w = 2.71 * discharge ** 0.557
    d = 0.349 * discharge ** 0.341
    r = d * w / (2 * d + w)
    v = 1 / n * r ** (2 / 3) * slope ** 0.5
    return v / 1000


def navigable_water_mask(dist):
    # Is this reach navigable by Kayak? From
    # [@rood2006instream,@zinke2018comparing] it seems that reaches with a
    # flow lower than 5m³/s are not navigable even by professional extreme
    # sport athletes, and generally even that number seems to be an outlier
    # with opinions starting at 8m³/s, so we take that as the cutoff.
    raster = rasterio.features.rasterize(
        (
            shp.shape
            for shp in RIVERS.shp.iterShapeRecords()
            if shp.record[3] > 0.9030899869919434
        ),  # log(8.)/log(10.)
        out_shape=dist.shape,
        transform=dist.transform,
    ).astype(bool)
    raster |= rasterio.features.rasterize(
        [OCEANS], out_shape=dist.shape, transform=dist.transform
    ).astype(bool)
    raster |= rasterio.features.rasterize(
        LAKES.iterShapes(), out_shape=dist.shape, transform=dist.transform
    ).astype(bool)
    raster |= rasterio.features.rasterize(
        (shp for shp in MORE_RIVERS.iterShapes() if shp.shapeTypeName != "NULL"),
        out_shape=dist.shape,
        transform=dist.transform,
    ).astype(bool)

    d_n, d_e, d_ne = [], [], []
    for y in range(1, dist.shape[0] + 1):
        (lon0, lat0) = dist.transform * (0, y)
        (lon1, lat1) = dist.transform * (1, y - 1)

        d = GEODESIC.inverse((lon0, lat0), [(lon0, lat1), (lon1, lat0), (lon1, lat1)])
        d_n.append(d[0, 0])
        d_e.append(d[1, 0])
        d_ne.append(d[2, 0])

    profile = dist.profile
    profile["dtype"] = rasterio.float64
    profile["count"] = 8
    fname = "adj{:}".format(dist.name)
    with rasterio.open(fname, "w", **profile) as dst:
        cell_distance = numpy.ones(dist.shape) * numpy.array(d_n)[:, None] / KAYAK_SPEED
        mask = raster[1:, :] & raster[:-1, :]
        # North
        with_rivers = dist.read(1)
        with_rivers[1:, :][mask] = cell_distance[1:, :][mask]
        dst.write(with_rivers.astype(rasterio.float64), 1)
        # South
        with_rivers = dist.read(5)
        with_rivers[:-1, :][mask] = cell_distance[:-1, :][mask]
        dst.write(with_rivers.astype(rasterio.float64), 5)

        cell_distance = numpy.ones(dist.shape) * numpy.array(d_e)[:, None] / KAYAK_SPEED
        mask = raster[:, 1:] & raster[:, :-1]
        # East
        with_rivers = dist.read(3)
        with_rivers[:, :-1][mask] = cell_distance[:, :-1][mask]
        dst.write(with_rivers.astype(rasterio.float64), 3)
        # West
        with_rivers = dist.read(7)
        with_rivers[:, 1:][mask] = cell_distance[:, 1:][mask]
        dst.write(with_rivers.astype(rasterio.float64), 7)

        cell_distance = (
            numpy.ones(dist.shape) * numpy.array(d_ne)[:, None] / KAYAK_SPEED
        )
        mask = raster[:-1, 1:] & raster[1:, :-1]
        # North-East
        with_rivers = dist.read(2)
        with_rivers[1:, :-1][mask] = cell_distance[1:, :-1][mask]
        dst.write(with_rivers.astype(rasterio.float64), 2)
        # South-West
        with_rivers = dist.read(6)
        with_rivers[:-1, 1:][mask] = cell_distance[:-1, 1:][mask]
        dst.write(with_rivers.astype(rasterio.float64), 6)

        mask = raster[1:, 1:] & raster[:-1, :-1]
        # South-East
        with_rivers = dist.read(4)
        with_rivers[:-1, :-1][mask] = cell_distance[:-1, :-1][mask]
        dst.write(with_rivers.astype(rasterio.float64), 4)
        # North-West
        with_rivers = dist.read(8)
        with_rivers[1:, 1:][mask] = cell_distance[1:, 1:][mask]
        dst.write(with_rivers.astype(rasterio.float64), 8)

    # At the end of the ``with rasterio.Env()`` block, context
    # manager exits and all drivers are de-registered.
    return rasterio.open(fname)


#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors:  Aric Hagberg <hagberg@lanl.gov>
#           Loïc Séguin-C. <loicseguin@gmail.com>
#           Dan Schult <dschult@colgate.edu>
#           Niels van Adrichem <n.l.m.vanadrichem@tudelft.nl>

from heapq import heappush, heappop
from itertools import count


def dijkstra_multisource(G, sources):
    """Uses Dijkstra's algorithm to find shortest weighted paths

    Parameters
    ----------
    G : [8, M, N] numpy array

    sources : non-empty iterable of nodes
        Starting nodes for paths. If this is just an iterable containing
        a single node, then all paths computed by this function will
        start from that node. If there are two or more nodes in this
        iterable, the computed paths may begin from any one of the start
        nodes.

    Returns
    -------
    distance : dictionary
        A mapping from node to shortest distance to that node from one
        of the source nodes.

    Raises
    ------
    NodeNotFound
        If any of `sources` is not in `G`.

    """
    push = heappush
    pop = heappop
    dist = numpy.full(G.shape[1:], numpy.inf, dtype=float)
    domain = numpy.full(G.shape[1:], -1, dtype=int)
    source_dist = {}
    seen = {}
    # fringe is heapq with 3-tuples (distance,c,node)
    # use the count c to avoid comparing nodes (may not be able to)
    c = count()
    fringe = []
    hexes = {}
    for i, (h, source) in enumerate(sources.items()):
        seen[source] = 0
        hexes[i] = h
        domain[source] = i
        push(fringe, (0, next(c), source))
    while fringe:
        print(len(fringe))
        (d, _, v) = pop(fringe)
        if numpy.isfinite(dist[v]):
            continue  # already searched this node.
        dist[v] = d
        for direction, e in zip(
            [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)],
            G[:, v[0], v[1]],
        ):
            if not numpy.isfinite(e):
                continue
            u = (v[0] + direction[0], v[1] + direction[1])
            if not 0 <= u[0] < G.shape[1]:
                continue
            if not 0 <= u[1] < G.shape[2]:
                continue
            vu_dist = dist[v] + e
            if u not in seen or vu_dist < seen[u]:
                if domain[u] != -1 and domain[u] != domain[v]:
                    i = min(domain[u], domain[v])
                    j = max(domain[u], domain[v])
                    source_dist[i, j] = min(
                        seen.get(u, numpy.inf) + vu_dist,
                        source_dist.get((i, j), numpy.inf),
                    )
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                domain[u] = domain[v]

    return hexes, domain, source_dist


AMERICAS = (
    -168.571541,  # w = xmin
    -56.028198,  # s = ymin
    -34.535395,  # e = xmax
    74.526716,  # n = ymax
)

BF = 500

for lon in range(-175, -45, 30):
    for lat in range(-60, 80, 20):
        print(lon, lat)

        fname = "{:}{:}.tif".format(lon, lat)
        try:
            ecoregions_file = ecoregion_tile_from_geocoordinates(lon, lat)
        except rasterio.errors.RasterioIOError:
            continue
        profile = ecoregions_file.profile
        height, width = ecoregions_file.shape
        profile["width"] = width + 2 * BF
        profile["height"] = height + 2 * BF
        profile["transform"] = ecoregions_file.transform
        profile["dtype"] = rasterio.float64
        profile["count"] = 8
        print(profile)

        try:
            distance_by_direction, areas = travel_time_raster(lon, lat, BF)
        except rasterio.errors.RasterioIOError:
            continue

        with rasterio.open(
            fname,
            "w",
            **profile,
        ) as dst:
            for i, band in enumerate(distance_by_direction.values(), 1):
                dst.write(band.astype(rasterio.float64), i)

        with open("areas-{:}.json".format(lat), "w") as areas_file:
            json.dump(areas, areas_file)

        travel_time = rasterio.open(fname)
        travel_time = navigable_water_mask(travel_time)
        try:
            travel_time = rasterio.open("adj{:}".format(fname))
        except rasterio.errors.RasterioIOError:
            continue

        land = sgeom.box(*travel_time.bounds).difference(OCEANS)
        xmin, ymin, xmax, ymax = land.bounds
        center = h3.geo_to_h3((ymin + ymax) / 2, (xmin + xmax) / 2, 5)
        corners = [h3.geo_to_h3(y, x, 5) for x in [xmin, xmax] for y in [ymin, ymax]]
        dist = 1
        for corner in corners:
            try:
                index = h3.h3_distance(center, corner)
            except SystemError:
                continue
            dist = max(dist, index + 1)
        pland = prep(land)
        dots = [
            hex
            for hex in h3.k_ring(center, dist)
            if pland.contains(sgeom.Point(h3.h3_to_geo(hex)[::-1]))
        ]

        back_transform = ~travel_time.transform

        def h3_to_rowcol(hex):
            lat, lon = h3.h3_to_geo(hex)
            col, row = back_transform * (lon, lat)
            return int(row), int(col)

        hexes, domain, dists = dijkstra_multisource(
            travel_time.read(), {dot: h3_to_rowcol(dot) for dot in dots}
        )
        fname = "adj{:}{:}-hexes.json".format(lon, lat)
        open(fname, "w").write(str(hexes))

        fname = "adj{:}{:}.json".format(lon, lat)
        json_dist = {}
        for (i, j), d in dists.items():
            i, j = hexes[i], hexes[j]
            i, j = min(i, j), max(i, j)
            json_dist.setdefault(int(i), {})[int(j)] = d
        with open(fname, "w") as d:
            json.dump(json_dist, d)

        profile = travel_time.profile
        profile["dtype"] = rasterio.int32
        profile["count"] = 1
        profile["count"] = 1
        fname = "vor{:}{:}.tif".format(lon, lat)
        with rasterio.open(fname, "w", **profile) as dst:
            dst.write(domain.astype(rasterio.int32), 1)

        ecoregions = ecoregions_file.read(1)
        assert domain[BF:-BF, BF:-BF].shape == ecoregions.shape
        c: t.DefaultDict[str, t.Counter[int]] = collections.defaultdict(
            collections.Counter
        )
        for voronoi_center, ecoregion, area in zip(
            domain[BF:-BF, BF:-BF].flat,
            ecoregions.flat,
            (numpy.array(areas[BF:-BF])[:, None] * numpy.ones(ecoregions.shape)).flat,
        ):
            c[int(hexes[voronoi_center])][int(ecoregion)] += area

        fname = "adj{:}{:}-regions.json".format(lon, lat)
        json.dump(c, open(fname, "w"))
