#!/home/cluster/gkaipi/.pyenv/shims/python
import json
import itertools
from pathlib import Path
import typing as t
import sys
import pickle
from itertools import count
from heapq import heappush as push, heappop as pop

import numpy
from tqdm import tqdm
from sqlalchemy import func
from sqlalchemy import select
from sqlalchemy.dialects.sqlite import insert

import shapely
import shapely.geometry as sgeom
from h3.api import basic_int as h3
from matplotlib import pyplot as plt

import rasterio
from raster_data import (
    tile_from_geocoordinates,
    gmted_tile,
    ecoregion_tile,
    boundingbox_from_tile,
    Tile,
)
from database import db
from ecoregions import TC
from earth import LAND, PLAND, BBOX, GEODESIC

DATABASE, TABLES = db()
RESOLUTION = 5

# ============
# Type aliases
# ============

RowCol = t.Tuple[int, int]
LonLat = t.Tuple[float, float]
H3Index = int

# =================================
# General geometry helper functions
# =================================
def as_point(h3_index: H3Index) -> sgeom.Point:
    y, x = h3.h3_to_geo(h3_index)
    return sgeom.Point(x, y)


def find_coast_hexagons() -> t.Set[H3Index]:
    coastal = set()
    # Coast resolution is 10m and coasts are not known to be straight, so we
    # can expect that every coastal hexagon contains at least one of the
    # coastline polygon coordinates.
    for geom in LAND.boundary.geoms:
        for x, y in geom.coords:
            coastal.add(h3.geo_to_h3(y, x, 5))
    return coastal


try:
    COAST = set(json.load(open("COAST.json")))
except (FileNotFoundError, json.JSONDecodeError):
    COAST = find_coast_hexagons()
    json.dump(list(COAST), open("COAST.json", "w"))


# =========================
# Travelling time functions
# =========================
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


# =================
# Raster operations
# =================
def prep_tile(tile: Tile) -> t.Tuple[numpy.ndarray, numpy.ndarray, rasterio.Affine]:
    """Load ecology data.

    Load the ecology metadata (elevations and ecoregions, in that order) of the
    tile (30°×20°) surrounding `containing`, plus a 500 pixel overlapping
    boundary region on each side, giving a (8200×5800) numpy array for the
    elevations (in m) and the ecoregions (integer categories)

    Returns
    =======
    Elevation: numpy.array
    Ecoregions: numpy.array
    transform: rasterio.Affine

    """
    print(f"Preparing rasters around for tile {tile}:")
    elevation_file = gmted_tile(tile)
    m_trafo = elevation_file.transform
    height, width = elevation_file.shape
    elevation = numpy.full((height + 1000, width + 1000), -100, int)
    ecoregions = numpy.full((height + 1000, width + 1000), 999, int)
    elevation[500:-500, 500:-500] = elevation_file.read(1)
    ecoregions[500:-500, 500:-500] = ecoregion_tile(tile).read(1)
    print("Loading margins of adjacent tiles…")
    west, south, east, north = boundingbox_from_tile(tile)
    try:
        nw = tile_from_geocoordinates(west - 15, north + 10)
        elevation[:500, :500] = gmted_tile(nw).read(1)[-500:, -500:]
        ecoregions[:500, :500] = ecoregion_tile(nw).read(1)[-500:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        n = tile_from_geocoordinates(west + 15, north + 10)
        elevation[:500, 500:-500] = gmted_tile(n).read(1)[-500:, :]
        ecoregions[:500, 500:-500] = ecoregion_tile(n).read(1)[-500:, :]
    except rasterio.RasterioIOError:
        pass
    try:
        ne = tile_from_geocoordinates(east + 15, north + 10)
        elevation[:500, -500:] = gmted_tile(ne).read(1)[-500:, :500]
        ecoregions[:500, -500:] = ecoregion_tile(ne).read(1)[-500:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        w = tile_from_geocoordinates(west - 15, south + 10)
        elevation[500:-500, :500] = gmted_tile(w).read(1)[:, -500:]
        ecoregions[500:-500, :500] = ecoregion_tile(w).read(1)[:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        e = tile_from_geocoordinates(east + 15, south + 10)
        elevation[500:-500, -500:] = gmted_tile(e).read(1)[:, :500]
        ecoregions[500:-500, -500:] = ecoregion_tile(e).read(1)[:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        sw = tile_from_geocoordinates(west - 15, south - 10)
        elevation[-500:, :500] = gmted_tile(sw).read(1)[:500, -500:]
        ecoregions[-500:, :500] = ecoregion_tile(sw).read(1)[:500, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        s = tile_from_geocoordinates(west + 15, south - 10)
        elevation[-500:, 500:-500] = gmted_tile(s).read(1)[:500, :]
        ecoregions[-500:, 500:-500] = ecoregion_tile(s).read(1)[:500, :]
    except rasterio.RasterioIOError:
        pass
    try:
        se = tile_from_geocoordinates(east + 15, south - 10)
        elevation[-500:, -500:] = gmted_tile(se).read(1)[:500, :500]
        ecoregions[-500:, -500:] = ecoregion_tile(se).read(1)[:500, :500]
    except rasterio.RasterioIOError:
        pass

    transform = rasterio.Affine(
        m_trafo.a,
        0,
        m_trafo.c - 500 * m_trafo.a,
        0,
        m_trafo.e,
        m_trafo.f - 500 * m_trafo.e,
    )
    print("Data loaded")
    return elevation, ecoregions, transform


def all_moore_neighbor_distances(
    elevation: numpy.array,
    transform: rasterio.Affine,
    terrain_coefficients: numpy.array,
) -> t.Dict[RowCol, numpy.array]:
    """Calculate the arrays of distances in all 8 neighbor directions.

    From an elevation raster and a terrain coefficient raster (higher
    coefficient = higher walking speed), localized using an affine
    transformation, compute the walking times 1 pixel into each direction (N,
    NE, E, …, NW).

    Return the resulting arrays inside a dictionary with index manipulations:

        {
         (-1, 0): distance to north,
         (-1, 1): distance to north east,
         (0, 1): distance to east,
         ...
         (-1, -1): distance to north west,
        }

    """
    # Compute the geodesic distances. They are constant for each row, which
    # corresponds to a constant latitude.
    d_n, d_e, d_ne = [], [], []
    for y in range(1, len(elevation) + 1):
        (lon0, lat0) = transform * (0, y)
        (lon1, lat1) = transform * (1, y - 1)

        d = GEODESIC.inverse((lon0, lat0), [(lon0, lat1), (lon1, lat0), (lon1, lat1)])
        d_n.append(d[0, 0])
        d_e.append(d[1, 0])
        d_ne.append(d[2, 0])
    distance_to_north = numpy.array(d_n)[:-1]
    slope_to_north = (
        100 * (elevation[1:, :] - elevation[:-1, :]) / distance_to_north[:, None]
    )
    tc_to_north = (terrain_coefficients[1:, :] + terrain_coefficients[:-1, :]) / 2
    north = distance_to_north[:, None] / (
        navigation_speed(slope_to_north) * tc_to_north
    )
    south = distance_to_north[:, None] / (
        navigation_speed(-slope_to_north) * tc_to_north
    )
    del distance_to_north, slope_to_north, tc_to_north

    distance_to_east = numpy.array(d_e)
    slope_to_east = (
        100 * (elevation[:, 1:] - elevation[:, :-1]) / distance_to_east[:, None]
    )
    tc_to_east = (terrain_coefficients[:, 1:] + terrain_coefficients[:, :-1]) / 2
    east = distance_to_east[:, None] / (navigation_speed(slope_to_east) * tc_to_east)
    west = distance_to_east[:, None] / (navigation_speed(-slope_to_east) * tc_to_east)
    del distance_to_east, slope_to_east, tc_to_east

    distance_to_northeast = numpy.array(d_ne)[:-1]
    slope_to_northeast = (
        100 * (elevation[1:, 1:] - elevation[:-1, :-1]) / distance_to_northeast[:, None]
    )
    tc_to_northeast = (
        terrain_coefficients[1:, 1:] + terrain_coefficients[:-1, :-1]
    ) / 2
    northeast = distance_to_northeast[:, None] / (
        navigation_speed(slope_to_northeast) * tc_to_northeast
    )
    southwest = distance_to_northeast[:, None] / (
        navigation_speed(-slope_to_northeast) * tc_to_northeast
    )
    del distance_to_northeast, slope_to_northeast, tc_to_northeast
    distance_to_northwest = numpy.array(d_ne)[:-1]
    slope_to_northwest = (
        100 * (elevation[1:, :-1] - elevation[:-1, 1:]) / distance_to_northwest[:, None]
    )

    tc_to_northwest = (
        terrain_coefficients[1:, :-1] + terrain_coefficients[:-1, 1:]
    ) / 2

    southeast = distance_to_northwest[:, None] / (
        navigation_speed(-slope_to_northwest) * tc_to_northwest
    )
    northwest = distance_to_northwest[:, None] / (
        navigation_speed(slope_to_northwest) * tc_to_northwest
    )
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
    }


def distances_and_cache(tile: Tile):
    fname = "distances-{:s}{:d}{:s}{:d}.tif".format(*tile)

    try:
        distance_raster = rasterio.open(fname)
        return {
            (-1, 0): distance_raster.read(1),
            (-1, 1): distance_raster.read(2),
            (0, 1): distance_raster.read(3),
            (1, 1): distance_raster.read(4),
            (1, 0): distance_raster.read(5),
            (1, -1): distance_raster.read(6),
            (0, -1): distance_raster.read(7),
            (-1, -1): distance_raster.read(8),
        }, distance_raster.transform
    except rasterio.errors.RasterioIOError:
        print("Cache miss")

    elevation, ecoregions, transform = prep_tile(tile)

    terrain_coefficient_raster = TC[ecoregions]
    distance_by_direction = all_moore_neighbor_distances(
        elevation, transform, terrain_coefficient_raster
    )
    profile = rasterio.profiles.DefaultGTiffProfile()
    profile["height"] = elevation.shape[0]
    profile["width"] = elevation.shape[1]
    profile["transform"] = transform
    profile["dtype"] = rasterio.float64
    profile["count"] = 8
    del elevation, ecoregions, terrain_coefficient_raster

    with rasterio.open(
        fname,
        "w",
        **profile,
    ) as dst:
        for i, band in enumerate(distance_by_direction.values(), 1):
            dst.write(band.astype(rasterio.float64), i)
    return distance_by_direction, transform


def distances_from_focus(
    source: t.Tuple[int, int],
    destinations: t.Optional[t.Set[RowCol]],
    distance_by_direction: t.Dict[RowCol, numpy.array],
    pred: t.Optional = None,
) -> numpy.array:
    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    d = distance_by_direction[0, 1]
    dist: numpy.array = numpy.full((d.shape[0], d.shape[1] + 1), numpy.inf, dtype=float)
    seen: t.Dict[t.Tuple[int, int], float] = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples
    fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []
    push(fringe, (0, next(c), source))

    def moore_neighbors(
        r0: int, c0: int
    ) -> t.Iterable[t.Tuple[t.Tuple[int, int], float]]:
        for (r, c), d in distance_by_direction.items():
            r1, c1 = r0 + r, c0 + c
            r = min(r0, r1)
            c = min(c0, c1)
            if 0 <= r < d.shape[0] and 0 <= c < d.shape[1]:
                yield (r1, c1), d[r, c]

    while fringe:
        (d, _, spot) = pop(fringe)
        if numpy.isfinite(dist[spot]):
            continue  # already searched this node.
        dist[spot] = d
        if destinations is not None and spot in destinations:
            destinations.remove(spot)
            if not destinations:
                break

        for u, cost in moore_neighbors(*spot):
            vu_dist = dist[spot] + cost
            if numpy.isfinite(dist[u]):
                if vu_dist < dist[u]:
                    raise ValueError(
                        "Contradictory paths found. Do you have negative weights?"
                    )
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if pred is not None:
                    pred[u] = spot
    return dist


def distances_trimmed(source, points, distance_by_direction, transform):
    rmin = min(r for r, c in points) - 2
    rmax = max(r for r, c in points) + 3
    cmin = min(c for r, c in points) - 2
    cmax = max(c for r, c in points) + 3
    points = [(r - rmin, c - cmin) for r, c in points]

    r0, c0 = source[0] - rmin, source[1] - cmin

    dist = {
        (n, e): d[
            rmin - min(n, 0) : rmax - max(0, n), cmin - min(e, 0) : cmax - max(0, e)
        ]
        for (n, e), d in distance_by_direction.items()
    }

    return rmin, cmin, distances_from_focus((r0, c0), points, dist, pred=None)


# ================
# Process hexagons
# ================
def find_land_hexagons():
    hexagons = set()
    for poly in LAND.geoms:
        d = shapely.geometry.mapping(poly)
        hexagons |= h3.polyfill_geojson(d, 5)
    return hexagons | COAST


try:
    ALL = set(json.load(open("ALL.json")))
except (FileNotFoundError, json.JSONDecodeError):
    ALL = find_land_hexagons()
    json.dump(list(ALL), open("ALL.json", "w"))


def core_point(hexbin, distance_by_direction, transform):
    def rowcol(latlon):
        lat, lon = latlon
        if lon > 170:
            # FIXME: We can and need to do this because we are working on the
            # Americas and the Americas only. The generic solution is more
            # difficult.
            lon = lon - 360
        col, row = ~transform * (lon, lat)
        return int(row), int(col)

    points = [rowcol(latlon) for latlon in h3.h3_to_geo_boundary(hexbin)]
    rmin = min(r for r, c in points) - 1
    rmax = max(r for r, c in points) + 1
    cmin = min(c for r, c in points) - 1
    cmax = max(c for r, c in points) + 1
    points = [(r - rmin, c - cmin) for r, c in points]

    dist = {
        (n, e): d[
            rmin - min(n, 0) : rmax - max(0, n), cmin - min(e, 0) : cmax - max(0, e)
        ]
        for (n, e), d in distance_by_direction.items()
    }
    dist[0, -1] = dist[0, 1] = dist[0, 1] + dist[0, -1]
    dist[-1, -1] = dist[1, 1] = dist[1, 1] + dist[-1, -1]
    dist[1, -1] = dist[-1, 1] = dist[-1, 1] + dist[1, -1]
    dist[-1, 0] = dist[1, 0] = dist[1, 0] + dist[-1, 0]

    border = points[:]
    for i in range(-1, len(points) - 1):
        r0, c0 = points[i]
        r1, c1 = points[i + 1]
        border.append(
            (
                round(0.5 * (r0 + r1)),
                round(0.5 * (c0 + c1)),
            )
        )

    all_dist = 0
    for r0, c0 in border:
        all_dist = all_dist + distances_from_focus((r0, c0), set(border), dist)

    (r0, c0) = numpy.unravel_index(numpy.argmin(all_dist), all_dist.shape)
    return (r0 + rmin, c0 + cmin)


def tile_core_points(tile: Tile, skip_existing=True, n=0):
    west, south, east, north = boundingbox_from_tile(tile)

    def tile_contains(hexagon: H3Index):
        lat, lon = h3.h3_to_geo(hexagon)
        return (west < lon < east) and (south < lat < north)

    hexes_by_tile = {hex for hex in ALL if tile_contains(hex)}
    if skip_existing:
        hexes_by_tile -= {
            n
            for n, in DATABASE.execute(
                select(TABLES["nodes"].c.node_id).where(
                    TABLES["nodes"].c.h3longitude >= west,
                    TABLES["nodes"].c.h3latitude >= south,
                    TABLES["nodes"].c.h3longitude < east,
                    TABLES["nodes"].c.h3latitude < north,
                )
            )
        }

    distance_by_direction, transform = distances_and_cache(tile)

    try:
        values = []
        for hexagon in tqdm(hexes_by_tile):
            hlat, hlon = h3.h3_to_geo(hexagon)
            row, col = core_point(hexagon, distance_by_direction, transform)
            lon, lat = transform * (col, row)
            n += 1
            values.append(
                {
                    "node_id": hexagon,
                    "longitude": lon,
                    "latitude": lat,
                    "h3longitude": hlon,
                    "h3latitude": hlat,
                    "coastal": (hexagon in COAST),
                    "short": n,
                }
            )
    finally:
        if values:
            DATABASE.execute(insert(TABLES["nodes"]).values(values).on_conflict_do_nothing())
    return n


def voronoi_and_neighbor_distances(tile, skip_existing=True):
    west, south, east, north = boundingbox_from_tile(tile)
    nodes = [
        (node, short, (lon, lat))
        for node, short, lon, lat in DATABASE.execute(
            select(
                TABLES["nodes"].c.node_id,
                TABLES["nodes"].c.short,
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
            ).where(
                TABLES["nodes"].c.longitude >= west,
                TABLES["nodes"].c.latitude >= south,
                TABLES["nodes"].c.longitude < east,
                TABLES["nodes"].c.latitude < north,
            )
        )
    ]

    distance_by_direction, transform = distances_and_cache(tile)
    fname_v = "voronoi-{:s}{:d}{:s}{:d}.tif".format(*tile)
    fname_d = "min_distances-{:s}{:d}{:s}{:d}.tif".format(*tile)
    try:
        voronoi_allocation = rasterio.open(fname_v).read(1)
        min_distances = rasterio.open(fname_d).read(1)
    except rasterio.errors.RasterioIOError:
        height, width = distance_by_direction[1, 1].shape
        voronoi_allocation = numpy.zeros((height + 1, width + 1), dtype=numpy.uint32)
        min_distances = numpy.full(
            (height + 1, width + 1), numpy.inf, dtype=numpy.float32
        )

    def rowcol(lonlat: LonLat):
        nlon, nlat = lonlat
        if nlon > 170:
            # FIXME: We can and need to do this because we are working on the
            # Americas and the Americas only. The generic solution is more
            # difficult.
            nlon = nlon - 360
        ncol, nrow = ~transform * (nlon, nlat)
        return (int(nrow), int(ncol))

    try:
        for node, short, (lon, lat) in tqdm(nodes):
            if node > 100000000:
                # Make sure we have a small int to refer to this node!
                if short is None:
                    raise ValueError("Short ID of a hex node must not be None")
                # Node is an h3 index, not a river reach index
                environment = h3.k_ring(node, 2)
            else:
                # Node is a river reach end point
                environment = h3.k_ring(h3.geo_to_h3(lat, lon, 5), 1)
                environment.add(node)

            extremities = DATABASE.execute(
                select(
                    TABLES["nodes"].c.longitude,
                    TABLES["nodes"].c.latitude,
                    TABLES["nodes"].c.h3longitude,
                    TABLES["nodes"].c.h3latitude,
                ).where(
                    TABLES["nodes"].c.longitude != None,
                    TABLES["nodes"].c.latitude != None,
                    TABLES["nodes"].c.node_id.in_(environment),
                )
            ).fetchall()

            x, y, hx, hy = zip(*extremities)

            source = rowcol((lon, lat))
            nnodes = []
            nrowcol = []
            neighbors = DATABASE.execute(
                select(
                    TABLES["nodes"].c.node_id,
                    TABLES["nodes"].c.longitude,
                    TABLES["nodes"].c.latitude,
                ).where(
                    min(x) <= TABLES["nodes"].c.longitude,
                    TABLES["nodes"].c.longitude <= max(x),
                    min(y) <= TABLES["nodes"].c.latitude,
                    TABLES["nodes"].c.latitude <= max(y),
                )
            ).fetchall() + [
                # Add the hexagon neighbors as virtual targets, to have a
                # uniform breadth
                (None, hlon, hlat)
                for hlon, hlat in zip(hx, hy)
                if hlon is not None and hlat is not None
            ]

            if skip_existing:
                already_known = {
                    n
                    for n, in DATABASE.execute(
                        select(TABLES["edges"].c.node2).where(
                            TABLES["edges"].c.node1 == node,
                            TABLES["edges"].c.source == "grid",
                        )
                    )
                }
                to_be_known = {n for n, _, _ in neighbors} - {None}
                if already_known >= to_be_known and (
                    (node < 100000000)
                    or (
                        min_distances[source] == 0.0
                        and voronoi_allocation[source] == short
                    )
                ):
                    continue
                print(already_known - to_be_known, to_be_known - already_known)

            for nnode, nlon, nlat in neighbors:
                nrowcol.append(rowcol((nlon, nlat)))
                nnodes.append(nnode)

            rmin, cmin, array = distances_trimmed(
                source, nrowcol, distance_by_direction, transform
            )

            values = [
                {
                    "node1": node,
                    "node2": nnode,
                    "source": "grid",
                    "travel_time": array[nr - rmin, nc - cmin],
                }
                for nnode, (nr, nc) in zip(nnodes, nrowcol)
                if nnode is not None
            ]
            DATABASE.execute(
                insert(TABLES["edges"]).values(values).on_conflict_do_nothing()
            )

            if node > 100000000:
                # Node is an h3 index, not a river reach index, so update the voronoi shapes around it
                rows, cols = array.shape
                voronoi_allocation[rmin : rmin + rows, cmin : cmin + cols][
                    min_distances[rmin : rmin + rows, cmin : cmin + cols] >= array
                ] = short
                min_distances[rmin : rmin + rows, cmin : cmin + cols] = numpy.fmin(
                    min_distances[rmin : rmin + rows, cmin : cmin + cols], array
                )
    finally:
        profile = rasterio.profiles.DefaultGTiffProfile()
        profile["height"] = voronoi_allocation.shape[0]
        profile["width"] = voronoi_allocation.shape[1]
        profile["transform"] = transform
        del profile["dtype"]
        profile["count"] = 1

        with rasterio.open(
            fname_v,
            "w",
            dtype=numpy.uint32,
            **profile,
        ) as dst:
            dst.write(voronoi_allocation, 1)
        with rasterio.open(
            fname_d,
            "w",
            dtype=numpy.float32,
            **profile,
        ) as dst:
            dst.write(min_distances, 1)


def measure_ecoregions(tile: Tile):
    fname_v = "voronoi-{:s}{:d}{:s}{:d}.tif".format(*tile)
    # Only read the non-overlapping part
    voronoi_allocation = rasterio.open(fname_v).read(
        1, window=((500, 5800 - 500), (500, 8200 - 500))
    )
    ecoregion_data = ecoregion_tile(tile)
    ecoregions = ecoregion_data.read(1)
    transform = ecoregion_data.transform
    areas = numpy.zeros(len(ecoregions))
    for y, _ in tqdm(enumerate(areas)):
        # Approximate grid cells by rectangles.
        (lon0, lat0) = transform * (0, y)
        (lon1, lat1) = transform * (1, y + 1)

        d = GEODESIC.inverse((lon0, lat0), [(lon0, lat1), (lon1, lat0)])
        areas[y] = d[0, 0] * d[0, 1]

    values = t.DefaultDict(t.Counter)
    for area, voronoi, eco in zip(areas, voronoi_allocation, ecoregions):
        for v, e in zip(voronoi, eco):
            if e != 999 and v != 0:
                values[int(v)][int(e)] += area
    print(values)
    return values


if __name__ == "__main__":
    if "rivers" in sys.argv:
        # Separate module, which imports this one! It creates the distances
        # caches, too.
        ...

    if "core" in sys.argv:
        n, = DATABASE.execute(select(func.max(TABLES["nodes"].c.short))).fetchone()
        if n is None:
            n = 0
        else:
            n = (n // 140000 + 1) * 140000
        for tile in itertools.product(
            ["N", "S"],
            [10, 30, 50, 70],
            ["E", "W"],
            [0, 30, 60, 90, 120, 150, 180],
        ):
            try:
                n = tile_core_points(tile, n=n)
            except rasterio.errors.RasterioIOError:
                print("Tile {:s}{:d}{:s}{:d} not found.".format(*tile))

    if "popdense" in sys.argv:
        # A modification of Tallavaara's R Markdown script.
        ...

    if "sea" in sys.argv:
        # Separate module. This can run in parallel to the following steps, it
        # only needs all core points to be defined, but does not interfere with
        # voronoi cells or grid distances.
        ...

    if "voronoi" in sys.argv or "distances" in sys.argv:
        for tile in itertools.product(
            ["N", "S"], [10, 30, 50, 70], ["E", "W"], [0, 30, 60, 90, 120, 150, 180]
        ):
            print("Working on tile {:s}{:d}{:s}{:d}…".format(*tile))
            try:
                voronoi_and_neighbor_distances(tile)
            except rasterio.errors.RasterioIOError:
                print("Tile {:s}{:d}{:s}{:d} not found.".format(*tile))

    if "stitch" in sys.argv:
        # Separate module, and quite fast. You could run it on every Voronoi
        # tile whose neighbors have been computed.
        ...

    if "areas" in sys.argv:
        # The other steps can be run per tile (with some artefacts on the tile
        # boundaries, if the neighboring tiles have not finished the previous
        # step), but this one really needs all voronoi computations to have
        # finished.
        values = t.DefaultDict(t.Counter)
        node_from_short = dict(list(
            DATABASE.execute(
                select(
                    TABLES["nodes"].c.short,
                    TABLES["nodes"].c.node_id,
                ).where(TABLES["nodes"].c.short != None)
            )
        ))
        for tile in itertools.product(
            ["N", "S"], [10, 30, 50, 70], ["E", "W"], [0, 30, 60, 90, 120, 150, 180]
        ):
            print("Working on tile {:s}{:d}{:s}{:d}…".format(*tile))
            try:
                for short, counter in measure_ecoregions(tile).items():
                    values[node_from_short.get(short, -short)].update(counter)
            except rasterio.errors.RasterioIOError:
                print("Tile {:s}{:d}{:s}{:d} not found.".format(*tile))
            json.dump(values, open("areas.json", "w"), indent=2)
        DATABASE.execute(
            insert(TABLES["ecology"]).values(
                [
                    {"node": node, "ecoregion": ecoregion, "area": area / 1000000}
                    for node, areas in values.items()
                    for ecoregion, area in areas.items()
                ]
            )
        )

    if "populations" in sys.argv:
        # Needs `popdense` and `areas`.
        for node, ecoregion, area, popdensity in DATABASE.execute(
            select(
                TABLES["ecology"].c.node,
                TABLES["ecology"].c.ecoregion,
                TABLES["ecology"].c.area,
                TABLES["nodes"].c.popdensity,
            ).select_from(TABLES["ecology"].join(TABLES["nodes"]))
        ):
            DATABASE.execute(
                TABLES["ecology"]
                .update()
                .values(population_capacity=popdensity * area)
                .where(
                    (TABLES["ecology"].c.node == node)
                    & (TABLES["ecology"].c.ecoregion == ecoregion)
                )
            )
TODO: change caching behaviour
TODO: change borders for distances
