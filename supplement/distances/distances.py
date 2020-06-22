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


from database import db
from ecoregions import ECOREGIONS


try:
    __file__
except NameError:
    __file__ = 'this'

# Define some constants
GEODESIC: geodesic.Geodesic = geodesic.Geodesic()
RESOLUTION: int = 5
AREA: float = h3.hex_area(RESOLUTION, "km^2")

# The resolution of our hexes should be about 450 km² following Gavin (2017)
assert h3.hex_area(RESOLUTION - 1, "km^2") > 450 > AREA

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


def tile_from_geocoordinates(
        lon: float, lat: float
) -> t.Tuple[t.Literal["N", "S"],
             int,
             t.Literal["E", "W"],
             int]:
    """Turn a Longitude/Latitude (i.e. x, y) pair into a tile index

    The index describes the South West corner of the tile, which has a width of
    30° (anchored at the 0° meridian) and a height of 20° (centered at the
    equator, so anchored at ±10°). This kind of index is used eg. for the GMTED
    tiles.

    >>> tile_from_geocoordinates(15, 40)
    ('N', 30, 'E', 0)
    >>> tile_from_geocoordinates(-75, -20)
    ('S', 30, 'W', 90)
    >>> tile_from_geocoordinates(-179, -17)
    ('S', 30, 'W', 180)

    """
    southwest_corner_lon = int(lon // 30) * 30
    southwest_corner_lat = int((lat - 10) // 20) * 20 + 10
    ew: t.Literal["E", "W"] = "W" if southwest_corner_lon < 0 else "E"
    ns: t.Literal["N", "S"] = "S" if southwest_corner_lat < 0 else "N"
    return ns, abs(southwest_corner_lat), ew, abs(southwest_corner_lon)


def gmted_tile_from_geocoordinates(lon: float, lat: float) -> str:
    """Path to the GMTED tile for given geocoordinates.

    The return type is explicitly a `str`, not a `pathlib.Path`, because it
    contains a GDAL virtual file system component.

    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010N30E000/
    >>> gmted_tile_from_geocoordinates(15, 40)[-78:]
    'elevation/GMTED2010/GMTED2010N30E000_150.zip/30n000e_20101117_gmted_med150.tif'

    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010S30W090/
    >>> gmted_tile_from_geocoordinates(-75, -20)[-78:]
    'elevation/GMTED2010/GMTED2010S30W090_150.zip/30s090w_20101117_gmted_med150.tif'

    """
    ns, lat, ew, lon = tile_from_geocoordinates(lon, lat)
    file = "{:02d}{:1s}{:03d}{:1s}".format(
        lat, ns.lower(), lon, ew.lower())
    tile = "GMTED2010{:1s}{:02d}{:1s}{:03d}".format(
        ns.upper(), lat, ew.upper(), lon)
    path = (Path(__file__).absolute().parent.parent /
            f"elevation/GMTED2010/{tile:}_150.zip")
    return f"/vsizip/{path:}/{file:}_20101117_gmted_med150.tif"


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


biomes = {eco: record[3] for eco, record in ECOREGIONS.records.items()}


terrain_coefficients = {
    'Boreal Forests/Taiga': 1.0,
    'Deserts & Xeric Shrublands': 1.0,
    'Flooded Grasslands & Savannas': 1.0,
    'Mangroves': 1.0,
    'Mediterranean Forests, Woodlands & Scrub': 1.0,
    'Montane Grasslands & Shrublands': 1.0,
    'N/A': 1.0,
    'Temperate Broadleaf & Mixed Forests': 1.0,
    'Temperate Conifer Forests': 1.0,
    'Temperate Grasslands, Savannas & Shrublands': 1.0,
    'Tropical & Subtropical Coniferous Forests': 1.0,
    'Tropical & Subtropical Dry Broadleaf Forests': 1.0,
    'Tropical & Subtropical Grasslands, Savannas & Shrublands': 1.0,
    'Tropical & Subtropical Moist Broadleaf Forests': 1.0,
    'Tundra': 1.0,
}


TC = {b: terrain_coefficients[biome]
      for b, biome in biomes.items()}


def run_on_one_tile(
        lon: float, lat: float,
        db: Engine,
        hex: sqlalchemy.Table, dist: sqlalchemy.Table, eco: sqlalchemy.Table
) -> None:
    ns, lat0, ew, lon0 = tile_from_geocoordinates(lon, lat)
    gdal_path = gmted_tile_from_geocoordinates(lon, lat)
    ecoregions_path = "../ecoregions/ECOREGIONS-{0:02d}{1}{2:03d}{3}_20101117_gmted_med150.tif".format(
        lat0, ns.lower(), lon0, ew.lower())
    try:
        with rasterio.open(gdal_path) as elevation:
            with rasterio.open(ecoregions_path) as ecology:
                return areas(elevation, ecology,
                             db, hex, dist, eco)
    finally:
        plt.show()


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


def areas(
        raster: rasterio.DatasetReader,
        ecoraster: rasterio.DatasetReader,
        db: sqlalchemy.Engine,
        hex: sqlalchemy.Table, dist: sqlalchemy.Table, eco: sqlalchemy.Table
) -> None:
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
    i, j = raster.shape
    lon, lat = raster.transform * (j/2, i/2)
    starts: t.List[Index] = [h3.geo_to_h3(lat, lon, RESOLUTION)]
    failed: t.Set[Index] = set()
    finished: t.Set[Index] = set()

    elevation = raster.read(1)
    ecoregions = ecoraster.read(1)

    assert raster.transform.b == 0
    assert raster.transform.d == 0

    assert numpy.allclose(ecoraster.transform.a, raster.transform.a)
    assert ecoraster.transform.b == 0
    assert numpy.allclose(ecoraster.transform.c, raster.transform.c)
    assert ecoraster.transform.d == 0
    assert numpy.allclose(ecoraster.transform.e, raster.transform.e)
    assert numpy.allclose(ecoraster.transform.f, raster.transform.f)

    def geo_dist(x0: int, y0: int, x1: int, y1: int, cache: t.Dict[t.Tuple[int, int, int], float]={}) -> float:
        if x0 > x1:
            x0, x1 = x1, x0

        if y0 > y1:
            y0, y1 = y1, y0

        try:
            return cache[x1-x0, y0, y1]
        except KeyError:
            (lon0, lat0) = raster.transform * (x0, y0)
            (lon1, lat1) = raster.transform * (x1, y1)

            d = GEODESIC.inverse((lon0, lat0), (lon1, lat1))[0, 0]
            cache[x1-x0, y0, y1] = d
            return d

    while starts:
        start = starts.pop(0)
        if start in failed or start in finished:
            continue
        lat, lon = h3.h3_to_geo(start)
        neighbors: t.List[Index] = list(h3.k_ring(start, 2))

        distance_known = dist.select([dist.c.hexbin2]).where(dist.c.hexbin1 == start)
        if not set(neighbors) - {k[0] for k in db.execute(distance_known).fetchall()}:
            finished.add(start)
            starts.extend(n for n in neighbors
                          if n not in finished
                          if n not in failed)
            find_most = eco.select([eco.c.ecoregion]).where(eco.c.hexbin == start).order_by(eco.c.frequency)
            most = db.scalar(find_most)
            print(start, most)
            plt.scatter(lon, lat, c=ECOREGIONS.records.get(most, {11: '0.9'})[11])
            continue

        try:
            db.execute(hex.insert({"hexbin": start,
                                   "longitude": lon,
                                   "latitude": lat}))
        except sqlalchemy.exc.IntegrityError:
            pass

        try:
            distances, eco_count = distances_from_focus(
                start,
                neighbors,
                elevation,
                geo_dist,
                biomes,
                raster.transform,
                ecoregions)

            n: float = sum(eco_count.values())
            e: t.Mapping[int, float] = {r: c/n for r, c in eco_count.items()}

            # For debugging: Output the results
            print(lon, lat, neighbors,
                  {ECOREGIONS.records.get(r, {1: None})[1]: c for r, c in e.items()})

            # For debugging: prepare a plot
            r: int
            r, _ = eco_count.most_common(1)[0]
            plt.scatter(lon, lat, c=ECOREGIONS.records.get(r, {11: '0.1'})[11])

            # Dump results into DB
            for id, freq in eco_count.items():
                try:
                    db.execute(eco.insert({"hexbin": start,
                                           "ecoregion": int(id),
                                           "frequency": freq}))
                except sqlalchemy.exc.IntegrityError:
                    pass

            for n, d in distances.items():
                exists = sqlalchemy.select([hex.c.hexbin]).where(hex.c.hexbin == n)
                if not db.execute(exists).fetchone():
                    lat, lon = h3.h3_to_geo(n)
                    db.execute(hex.insert({"hexbin": n,
                                           "longitude": lon,
                                           "latitude": lat}))
                    if n not in starts:
                        starts.append(n)
                try:
                    db.execute(dist.insert({
                        "hexbin1": start,
                        "hexbin2": n,
                        "flat_distance": geographical_distance(start, n),
                        "distance": d,
                        "source": 0}))
                except sqlalchemy.exc.IntegrityError:
                    pass
            finished.add(start)

        except (IndexError, NotImplementedError):
            print("failed")
            failed.add(start)
    print(failed)


def distances_from_focus(
        i: Index,
        neighbors: t.List[Index],
        data: numpy.array,
        geo_dist: t.Callable[[int, int, int, int], float],
        biomes: t.Dict[int, str],
        transform: rasterio.Affine,
        ecoregions: rasterio.DatasetReader,
) -> t.Tuple[t.Dict[Index, float],
             t.Counter[int]]:
    origin = neighbors.index(i)
    points = numpy.array([h3.h3_to_geo(n)[::-1] for n in neighbors])
    indices = numpy.round(~transform * points.T).astype(int).T
    source = tuple(indices[origin])
    destinations = {tuple(j) for j in indices}

    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    dist: t.Dict[t.Tuple[int, int], float] = {}  # dictionary of final dists
    seen = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples
    fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []
    push(fringe, (0, next(c), source))

    e: t.Counter[int] = t.Counter()

    def moore_neighbors(
            x0: int, y0: int
    ) -> t.Iterable[t.Tuple[t.Tuple[int, int], float]]:
        elevation = data[y0, x0]
        # Look up terrain factor
        tf0 = TC.get(ecoregions[y0, x0], 1.0)
        for i, j in [(-1, -1), (-1, 0), (-1, 1),
                     (0, -1), (0, 1),
                     (1, -1), (1, 0), (1, 1)]:
            x1, y1 = x0 + i, y0 + j
            if x1 < 0 or y1 < 0 or x1 >= data.shape[1] or y1 >= data.shape[0]:
                raise IndexError("Attempted to move outside tile")
            horizontal = geo_dist(x0, y0, x1, y1)
            ed = data[y1, x1] - elevation
            tf1 = TC.get(ecoregions[y1, x1], 1.0)
            yield (x1, y1), horizontal / navigation_speed(
                ed / horizontal * 100) * (tf0 + tf1) / 2

    while fringe:
        (d, _, (x0, y0)) = pop(fringe)
        if (x0, y0) in dist:
            continue  # already searched this node.
        dist[(x0, y0)] = d
        if (x0, y0) in destinations:
            destinations.remove((x0, y0))
            if not destinations:
                break

        lon0, lat0 = transform * (x0, y0)
        j = h3.geo_to_h3(lat0, lon0, RESOLUTION)
        eco = ecoregions[y0, x0]
        if j == i:  # Looking at a tile within the focus hex
            e[eco] += numpy.cos(lat0 * numpy.pi / 180)

        if eco is None:
            print("No biome found at", lon0, lat0)
            continue

        for u, cost in moore_neighbors(x0, y0):
            vu_dist = dist[(x0, y0)] + cost
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
            elif vu_dist == seen[u]:
                pass
    return {n: dist[x, y]/3600 for n, (x, y) in zip(neighbors, indices)}, e


def plot_distances(db: sqlalchemy.engine.Connectable, dist: sqlalchemy.Table) -> None:
    distances = dist.select([dist.c.flat_distance, dist.c.distance])
    x, y = zip(*db.execute(distances))
    plt.scatter(x, y)


if __name__ == '__main__':
    # FIXME: Use Argparser instead
    import sys
    db, tables = db(sys.argv[3])
    run_on_one_tile(float(sys.argv[2]), float(sys.argv[1]), db, **tables)
    plot_distances(db, tables['dist'])
