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

# [@pozzi2008accessibility]: Walking speeds
# Open or sparse grasslands, croplands (> 50%, or with open woody vegetation, or irrigated), mosaic of forest/croplands or forest/savannah, urban areas: 3 km/h
# Deciduous shrubland or woodland, closed grasslands, tree crops, desert (sandy or stony) and dunes, bare rock: 1.5 km/h
# Lowland forest (deciduous or degraded evergreen), swamp bushland and grassland, salt hardpans: 1 km/h
# Submontane and montane forest: 0.6 km/h
# Closed evergreen lowland forest, swamp forest, mangrove: 0.3 km/h

# The off-road navigation speeds are for temperate forests, so they get a terrain factor of 1.
terrain_coefficients = {
    'Temperate Broadleaf & Mixed Forests': 1.0,
    'Temperate Conifer Forests': 1.0,
    'Deserts & Xeric Shrublands': 1.5,
    'Boreal Forests/Taiga': 1.0,
    'Flooded Grasslands & Savannas': 1.0,
    'Mangroves': 0.3,
    'Mediterranean Forests, Woodlands & Scrub': 1.0,
    'Montane Grasslands & Shrublands': 1.5,
    'N/A': 0.05,
    'Temperate Grasslands, Savannas & Shrublands': 3.0,
    'Tropical & Subtropical Coniferous Forests': 1.0,
    'Tropical & Subtropical Dry Broadleaf Forests': 1.0,
    'Tropical & Subtropical Grasslands, Savannas & Shrublands': 3.0,
    'Tropical & Subtropical Moist Broadleaf Forests': 1.0,
    'Tundra': 1.5,
}


TC = numpy.array([
    terrain_coefficients[biomes.get(b, "N/A")]
    for b in range(1000)])


def run_on_one_tile(
        lon: float, lat: float,
        db: sqlalchemy.engine.Engine,
        hex: sqlalchemy.Table, dist: sqlalchemy.Table, eco: sqlalchemy.Table
) -> None:
    ns, lat0, ew, lon0 = tile_from_geocoordinates(lon, lat)
    gdal_path = gmted_tile_from_geocoordinates(lon, lat)
    ecoregions_path = "../ecoregions/ECOREGIONS-{0:02d}{1}{2:03d}{3}_20101117_gmted_med150.tif".format(
        lat0, ns.lower(), lon0, ew.lower())
    with rasterio.open(gdal_path) as elevation:
        with rasterio.open(ecoregions_path) as ecology:
            return areas(elevation, ecology,
                         db, hex, dist, eco)


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
) -> numpy.array:
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
    distance_to_east = numpy.array(d_e)
    distance_to_northeast = numpy.array(d_ne)[:-1]
    distance_to_northwest = distance_to_northeast

    slope_to_north = 100 * (elevation[1:, :] - elevation[:-1, :]) / distance_to_north[:, None]
    slope_to_east = 100 * (elevation[:, 1:] - elevation[:, :-1]) / distance_to_east[:, None]
    slope_to_northeast = 100 * (elevation[1:, 1:] - elevation[:-1, :-1]) / distance_to_northeast[:, None]
    slope_to_northwest = 100 * (elevation[1:, :-1] - elevation[:-1, 1:]) / distance_to_northwest[:, None]

    north = distance_to_north[:, None] / navigation_speed(slope_to_north)
    northeast = distance_to_northeast[:, None] / navigation_speed(slope_to_northeast)
    east = distance_to_east[:, None] / navigation_speed(slope_to_east)
    southeast = distance_to_northwest[:, None] / navigation_speed(-slope_to_northwest)
    south = distance_to_north[:, None] / navigation_speed(-slope_to_north)
    southwest = distance_to_northeast[:, None] / navigation_speed(-slope_to_northeast)
    west = distance_to_east[:, None] / navigation_speed(-slope_to_east)
    northwest = distance_to_northwest[:, None] / navigation_speed(slope_to_northwest)

    return {
        (0, -1): north,
        (1, -1): northeast,
        (1, 0): east,
        (1, 1): southeast,
        (0, 1): south,
        (-1, 1): southwest,
        (-1, 0): west,
        (-1, -1): northwest
    }


def areas(
        raster: rasterio.DatasetReader,
        ecoraster: rasterio.DatasetReader,
        db: sqlalchemy.engine.Engine,
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

    for h, eco_count in hex_ecoregions(ecoregions, ecoraster.transform).items():
        lat, lon = h3.h3_to_geo(h)
        try:
            db.execute(hex.insert({"hexbin": h,
                                   "longitude": lon,
                                   "latitude": lat}))
        except sqlalchemy.exc.IntegrityError:
            pass
        for id, freq in eco_count.items():
            try:
                db.execute(eco.insert({"hexbin": h,
                                       "ecoregion": int(id),
                                       "frequency": freq}))
            except sqlalchemy.exc.IntegrityError:
                pass
    del lat, lon

    distance_by_direction = all_pairwise_distances(
        elevation, raster.transform)

    terrain_coefficient_raster = TC[ecoregions]

    while starts:
        start = starts.pop(0)
        if start in failed or start in finished:
            continue
        neighbors: t.List[Index] = list(h3.k_ring(start, 2))

        distance_known = sqlalchemy.select([dist.c.hexbin2]).where(dist.c.hexbin1 == start)
        if not set(neighbors) - {k[0] for k in db.execute(distance_known).fetchall()}:
            finished.add(start)
            starts.extend(n for n in neighbors
                          if n not in finished
                          if n not in failed)
            find_most = sqlalchemy.select([eco.c.ecoregion]).where(eco.c.hexbin == start).order_by(eco.c.frequency)
            most = db.scalar(find_most)
            print(start, most)
            continue

        try:
            origin = neighbors.index(start)
            points = numpy.array([h3.h3_to_geo(n)[::-1] for n in neighbors])
            indices = numpy.round(~raster.transform * points.T).astype(int).T
            destinations: t.List[t.Tuple[int, int]] = [
                (i, j) for (i, j) in indices]
            source: t.Tuple[int, int] = destinations[origin]

            distances = distances_from_focus(
                source,
                set(destinations),
                distance_by_direction,
                terrain_coefficient_raster)

            local_d = {n: distances[t[1], t[0]] for n, t in zip(neighbors, destinations)}
            # For debugging: Output the results
            print(local_d)

            for n, d in local_d.items():
                if n not in starts and n not in failed and n not in finished:
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


def hex_ecoregions(
        ecoregions: numpy.array,
        transform: rasterio.Affine
) -> t.Dict[h3.H3Index, t.Counter[int]]:
    c: t.Dict[h3.H3Index, t.Counter[int]] = t.DefaultDict(t.Counter)
    for y, row in enumerate(ecoregions):
        (_, lat) = transform * (0, y)
        area = numpy.cos(lat * numpy.pi / 180)
        for x, eco in enumerate(row):
            (lon, lat) = transform * (x, y)
            index: h3.H3Index = h3.geo_to_h3(lat, lon, RESOLUTION)
            c[index][eco] += area
    return c


def distances_from_focus(
        source: t.Tuple[int, int],
        destinations: t.Set[t.Tuple[int, int]],
        distance_by_direction: t.Dict[t.Tuple[int, int], numpy.array],
        terrain_coefficients_raster: numpy.array,
) -> t.Dict[t.Tuple[int, int], float]:

    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    dist: numpy.array = numpy.full(
        terrain_coefficients_raster.shape,
        numpy.nan, dtype=float)
    seen = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples
    fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []
    push(fringe, (0, next(c), source))

    def moore_neighbors(
            x0: int, y0: int
    ) -> t.Iterable[t.Tuple[t.Tuple[int, int], float]]:
        tf0 = terrain_coefficients_raster[y0, x0]
        for (i, j), d in distance_by_direction.items():
            x1, y1 = x0 + i, y0 + j
            if x1 < 0 or y1 < 0:
                raise IndexError("Attempted to move outside tile")
            tf1 = terrain_coefficients_raster[y1, x1]
            yield (x1, y1), d[y0, x0] / (tf0 + tf1) * 2

    while fringe:
        (d, _, (x0, y0)) = pop(fringe)
        if numpy.isfinite(dist[y0, x0]):
            continue  # already searched this node.
        dist[y0, x0] = d
        if (x0, y0) in destinations:
            destinations.remove((x0, y0))
            if not destinations:
                break

        for u, cost in moore_neighbors(x0, y0):
            vu_dist = dist[y0, x0] + cost
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
            elif vu_dist == seen[u]:
                pass
    return dist / 3600


def plot_distances(db: sqlalchemy.engine.Connectable, dist: sqlalchemy.Table) -> None:
    distances = sqlalchemy.select([dist.c.flat_distance, dist.c.distance])
    x, y = zip(*db.execute(distances))
    plt.scatter(x, y, marker='x', s=40, alpha=0.2)


if __name__ == '__main__':
    # FIXME: Use Argparser instead
    import sys
    engine, tables = db(sys.argv[3])
    t_hex = tables["hex"]
    t_dist = tables["dist"]
    t_eco = tables["eco"]
    run_on_one_tile(float(sys.argv[2]), float(sys.argv[1]),
                    engine, hex=t_hex, dist=t_dist, eco=t_eco)
    plot_distances(engine, t_dist)
