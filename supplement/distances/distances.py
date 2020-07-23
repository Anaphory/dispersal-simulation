"""
distances – various distance measures for prehistoric anthropology
"""

import json
import numpy
import zipfile
import collections
import typing as t
from pathlib import Path
from itertools import count
from dataclasses import dataclass
from heapq import heappush as push, heappop as pop

import numpy
import sqlalchemy
from sqlalchemy.sql import func

import matplotlib.pyplot as plt
import matplotlib.collections
import matplotlib.transforms as mtransforms

import cartopy.geodesic as geodesic
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

import shapefile
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

import rasterio

from h3 import h3
from h3.h3 import H3Index as Index

import networkx

from database import db
from ecoregions import ECOREGIONS, TC, RESOLUTION

try:
    __file__
except NameError:
    __file__ = 'this'

ArrayIndex = t.Tuple[int, int]

# Define some constants
GEODESIC: geodesic.Geodesic = geodesic.Geodesic()

P = t.TypeVar("P", bound=sgeom.Point)


class Point(sgeom.Point):  # type: ignore
    @classmethod
    def from_h3(cls: t.Type[P], index: Index) -> P:
        lat, lon = h3.h3_to_geo(index)
        return cls(lon, lat)

    def __repr__(self) -> str:
        return f"Point(longitude={self.x:}, latitude={self.y:})"

    @property
    def longitude(self) -> float:
        return self.x

    @property
    def latitude(self) -> float:
        return self.y


land_shp_fname = shpreader.natural_earth(
    resolution='50m', category='physical', name='land')

land_geom = unary_union(
    [record.geometry
     for record in shpreader.Reader(land_shp_fname).records()
     if record.attributes.get('featurecla') != "Null island"])

LAND = prep(land_geom)

AMERICAS = (
    -168.571541,  # w = xmin
    -56.028198,   # s = ymin
    -34.535395,   # e = xmax
    74.526716,    # n = ymax
)

from raster_data import tile_from_geocoordinates, ecoregion_tile_from_geocoordinates, gmted_tile_from_geocoordinates

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
    return 0.11 + 0.67 * numpy.exp(-(slope + 2.0) ** 2 / 1800.)


def ring_in_order(center: h3.H3Index, k: int) -> t.Iterable[h3.H3Index]:
    """Return a hexagonal approximation of a circle.

    Give the elements of the k-ring around center in clockwise or anticlockwise
    order.

    """
    points: t.Set[h3.H3Index] = h3.k_ring_distances(center, k)[k]
    start = points.pop()
    n = start
    yield start
    while points:
        n = (h3.k_ring_distances(n, 1)[1] & points).pop()
        points.remove(n)
        yield n

def run_on_one_tile(
        lon: float, lat: float,
        db: sqlalchemy.engine.Engine,
        hex: sqlalchemy.Table, t_dist: sqlalchemy.Table
) -> t.Set[ArrayIndex]:
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
    elevation = numpy.full((height + 1000, width + 1000), -100, int)
    ecoregions = numpy.full((height + 1000, width + 1000), 999, int)
    elevation[500:-500, 500:-500] = elevation_file.read(1)
    ecoregions[500:-500, 500:-500] = ecoregion_tile_from_geocoordinates(lon, lat).read(1)
    print("Loading adjacent data…")
    try:
        elevation[:500, :500] = (gmted_tile_from_geocoordinates(lon - 30, lat + 20)).read(1)[-500:, -500:]
        ecoregions[:500, :500] = ecoregion_tile_from_geocoordinates(lon - 30, lat + 20).read(1)[-500:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[:500, 500:-500] = (gmted_tile_from_geocoordinates(lon, lat + 20)).read(1)[-500:, :]
        ecoregions[:500, 500:-500] = ecoregion_tile_from_geocoordinates(lon, lat + 20).read(1)[-500:, :]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[:500, -500:] = (gmted_tile_from_geocoordinates(lon + 30, lat + 20)).read(1)[-500:, :500]
        ecoregions[:500, -500:] = ecoregion_tile_from_geocoordinates(lon + 30, lat + 20).read(1)[-500:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[500:-500, :500] = (gmted_tile_from_geocoordinates(lon - 30, lat)).read(1)[:, -500:]
        ecoregions[500:-500, :500] = ecoregion_tile_from_geocoordinates(lon - 30, lat).read(1)[:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[500:-500, -500:] = (gmted_tile_from_geocoordinates(lon + 30, lat)).read(1)[:, :500]
        ecoregions[500:-500, -500:] = ecoregion_tile_from_geocoordinates(lon + 30, lat).read(1)[:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, :500] = (gmted_tile_from_geocoordinates(lon - 30, lat - 20)).read(1)[:500, -500:]
        ecoregions[-500:, :500] = ecoregion_tile_from_geocoordinates(lon - 30, lat - 20).read(1)[:500, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, 500:-500] = (gmted_tile_from_geocoordinates(lon, lat - 20)).read(1)[:500, :]
        ecoregions[-500:, 500:-500] = ecoregion_tile_from_geocoordinates(lon, lat - 20).read(1)[:500, :]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, -500:] = (gmted_tile_from_geocoordinates(lon + 30, lat - 20)).read(1)[:500, :500]
        ecoregions[-500:, -500:] = ecoregion_tile_from_geocoordinates(lon + 30, lat - 20).read(1)[:500, :500]
    except rasterio.RasterioIOError:
        pass

    print("Computing hex extents…")
    transform = rasterio.Affine(m_trafo.a, 0, m_trafo.c - 500 * m_trafo.a,
                                0, m_trafo.e, m_trafo.f - 500 * m_trafo.e)
    def rowcol(latlon):
        lat, lon = latlon
        if lon > 170:
            # FIXME: We can and need to do this because we are working on the
            # Americas and the Americas only. The generic solution is more
            # difficult.
            lon = lon - 360
        col, row = ~transform * (lon, lat)
        return int(row), int(col)

    center: t.Dict[h3.H3Index, t.Tuple[int, int]] = {}
    cs = sqlalchemy.select([hex.c.hexbin, hex.c.vlongitude, hex.c.vlatitude]).where(hex.c.vlongitude != None)
    for h, lon, lat in db.execute(cs).fetchall():
        row, col = rowcol((lat, lon))
        if 500 <= col < width + 500 and 500 < row < height + 500:
            center[h] = (row, col)

    print("Computing terrain coefficients…")
    terrain_coefficient_raster = TC[ecoregions]

    print("Computing distances on the grid…")
    distance_by_direction = all_pairwise_distances(
        elevation, transform, terrain_coefficient_raster)

    failed: t.Set[ArrayIndex] = set()
    finished: t.Set[ArrayIndex] = set()

    for start, (row, col) in center.items():
        print(f"Computing area around {start:}…")
        lon, lat = transform * (col, row)
        hexbin = h3.geo_to_h3(lat, lon, RESOLUTION)
        this, neighbors1, neighbors2, neighbors3 = list(h3.k_ring_distances(hexbin, 3))
        neighbors: t.Set[h3.H3Index] = neighbors1 | neighbors2
        distance_known = sqlalchemy.select([t_dist.c.hexbin2]).where(
            (t_dist.c.hexbin1 == start) & (t_dist.c.source == 0))
        missing = (neighbors - {k[0] for k in db.execute(distance_known).fetchall()}) & set(center)
        if not missing:
            print("Already known.")
            finished.add(start)
            continue
        print("Not found:", missing)
        if start in failed or start in finished:
            continue

        points = []
        for n in neighbors2 | neighbors3:
            try:
                points.append(center[n])
            except KeyError:
                points.append(rowcol(h3.h3_to_geo(n)))
        rmin = min(p[0] for p in points)
        rmax = max(p[0] for p in points)
        cmin = min(p[1] for p in points)
        cmax = max(p[1] for p in points)
        assert rmin >= 0
        assert cmin >= 0

        def rel_rowcol(latlon):
            lat, lon = latlon
            if lon > 170:
                lon = lon - 360
            col, row = ~transform * (lon, lat)
            return int(row) - rmin, int(col) - cmin

        destinations = [(center[n][0] - rmin, center[n][1] - cmin)
                        for n in neighbors
                        if n in center]

        print(f"Looking for {neighbors:}…")
        distances = distances_from_focus(
            (center[start][0] - rmin, center[start][1] - cmin),
            set(destinations),
            {(n, e): d[rmin-min(n, 0):rmax + 1 - max(0, n),
                       cmin-min(e, 0):cmax + 1 - max(0, e)]
             for (n, e), d in distance_by_direction.items()},
            pred=None
        )

        # For debugging: Output the results
        distances_by_center = {}
        for t, d in center.items():
            try:
                dist = distances[d[0] - rmin, d[1] - cmin]
                if numpy.isfinite(dist):
                    distances_by_center[t] = dist
            except IndexError:
                continue

        with db.begin() as conn:
            for n, d in distances_by_center.items():
                try:
                    conn.execute(t_dist.insert({
                        "hexbin1": start,
                        "hexbin2": n,
                        "distance": d,
                        "source": 0}).prefix_with('OR IGNORE'))
                except sqlalchemy.exc.IntegrityError:
                    pass
        finished.add(start)

    return failed


def geographical_distance(index1: Index, index2: Index) -> float:
    """Calculate the geodesic distance between two hex centers, in meters.

    >>> i = h3.geo_to_h3(
    ...   numpy.random.random() * 180 - 90,
    ...   numpy.random.random() * 360 - 180,
    ...   RESOLUTION
    ...   )
    >>> i_, neighbors = h3.k_ring_distances(i, 1)
    >>> i_ == {i}
    True
    >>> edge = h3.edge_length(RESOLUTION, 'm')
    >>> ex = 3**0.5 * edge
    >>> ex
    14799.349254643997
    >>> for n in neighbors:
    ...   print(0.7*ex < geographical_distance(i, n) < 1.3*ex
    ...     or geographical_distance(i, n))
    True
    True
    True
    True
    True
    True
    """
    lat1, lon1 = h3.h3_to_geo(index1)
    lat2, lon2 = h3.h3_to_geo(index2)
    return numpy.asarray(GEODESIC.inverse(
        (lon1, lat1),
        (lon2, lat2)))[0, 0]


def all_pairwise_distances(
        elevation: numpy.array,
        transform: rasterio.Affine,
        terrain_coefficients: numpy.array,
) -> t.Dict[ArrayIndex, numpy.array]:
    d_n, d_e, d_ne = [], [], []
    for y in range(1, len(elevation) + 1):
        (lon0, lat0) = transform * (0, y)
        (lon1, lat1) = transform * (1, y - 1)

        d = GEODESIC.inverse((lon0, lat0), [
            (lon0, lat1), (lon1, lat0), (lon1, lat1)])
        d_n.append(d[0, 0])
        d_e.append(d[1, 0])
        d_ne.append(d[2, 0])
    distance_to_north = numpy.array(d_n)[:-1]
    slope_to_north = 100 * (elevation[1:, :] - elevation[:-1, :]) / distance_to_north[:, None]
    tc_to_north = (terrain_coefficients[1:, :] + terrain_coefficients[:-1, :]) / 2
    north = distance_to_north[:, None] / (navigation_speed(slope_to_north) * tc_to_north)
    south = distance_to_north[:, None] / (navigation_speed(-slope_to_north) * tc_to_north)
    del distance_to_north, slope_to_north, tc_to_north

    distance_to_east = numpy.array(d_e)
    slope_to_east = 100 * (elevation[:, 1:] - elevation[:, :-1]) / distance_to_east[:, None]
    tc_to_east = (terrain_coefficients[:, 1:] + terrain_coefficients[:, :-1]) / 2
    east = distance_to_east[:, None] / (navigation_speed(slope_to_east) * tc_to_east)
    west = distance_to_east[:, None] / (navigation_speed(-slope_to_east) * tc_to_east)
    del distance_to_east, slope_to_east, tc_to_east

    distance_to_northeast = numpy.array(d_ne)[:-1]
    slope_to_northeast = 100 * (elevation[1:, 1:] - elevation[:-1, :-1]) / distance_to_northeast[:, None]
    tc_to_northeast = (terrain_coefficients[1:, 1:] + terrain_coefficients[:-1, :-1]) / 2
    northeast = distance_to_northeast[:, None] / (navigation_speed(slope_to_northeast) * tc_to_northeast)
    southwest = distance_to_northeast[:, None] / (navigation_speed(-slope_to_northeast) * tc_to_northeast)
    del distance_to_northeast, slope_to_northeast, tc_to_northeast
    distance_to_northwest = numpy.array(d_ne)[:-1]
    slope_to_northwest = 100 * (elevation[1:, :-1] - elevation[:-1, 1:]) / distance_to_northwest[:, None]

    tc_to_northwest = (terrain_coefficients[1:, :-1] + terrain_coefficients[:-1, 1:]) / 2

    southeast = distance_to_northwest[:, None] / (navigation_speed(-slope_to_northwest) * tc_to_northwest)
    northwest = distance_to_northwest[:, None] / (navigation_speed(slope_to_northwest) * tc_to_northwest)
    del distance_to_northwest, slope_to_northwest, tc_to_northwest

    return {
        (-1, 0): north,
        (-1, 1): northeast,
        (0, 1): east,
        (1, 1): southeast,
        (1, 0): south,
        (1, -1): southwest,
        (0, -1): west,
        (-1, -1): northwest
    }



def distances_from_focus(
        source: ArrayIndex,
        destinations: t.Optional[t.Set[ArrayIndex]],
        distance_by_direction: t.Dict[ArrayIndex, numpy.array],
        pred: t.Optional[t.Dict[Index, Index]]
) -> numpy.array:
    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    d = distance_by_direction[0, 1]
    dist: numpy.array = numpy.full((d.shape[0], d.shape[1] + 1),
        numpy.nan, dtype=float)
    seen: t.Dict[ArrayIndex, float] = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples
    fringe: t.List[t.Tuple[float, int, ArrayIndex]] = []
    push(fringe, (0, next(c), source))

    def moore_neighbors(
            r0: int, c0: int
    ) -> t.Iterable[t.Tuple[ArrayIndex, float]]:
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
                    raise ValueError('Contradictory paths found:',
                                    'negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if pred is not None:
                    pred[u] = spot
    return dist


def plot_distances(db: sqlalchemy.engine.Connectable, dist: sqlalchemy.Table) -> None:
    distances = sqlalchemy.select([dist.c.flat_distance, dist.c.distance])
    x, y = zip(*db.execute(distances))
    plt.scatter(x, y, marker='x', s=40, alpha=0.2)


def lonlat(s: str):
    lon, lat = s.split(",")
    return float(lon), float(lat)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("--tile", action="append", default=[], type=lonlat)
    args = parser.parse_args()
    engine, tables = db(args.database)
    t_hex = tables["hex"]
    t_dist = tables["dist"]
    t_eco = tables["eco"]

    tiles = args.tile or [(lon, lat)
             for lon in range(-165, 195, 30)
             for lat in range(60, -90, -20)]
    for lon, lat in tiles:
        try:
            run_on_one_tile(lon, lat, engine, t_hex, t_dist)
        except rasterio.RasterioIOError:
            continue

