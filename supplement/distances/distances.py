"""
distances – various distance measures for prehistoric anthropology
"""

import json
import numpy
import zipfile
import collections
import typing as t
from pathlib import Path
from dataclasses import dataclass

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

class Point(sgeom.Point):
    @classmethod
    def from_h3(cls, index: h3.H3Index):
        lat, lon = h3.h3_to_geo(index)
        return cls(lon, lat)

    def __repr__(self):
        return f"Point(longitude={self.x:}, latitude={self.y:})"

    @property
    def longitude(self):
        return self.x

    @property
    def latitude(self):
        return self.y

# Define classes
@dataclass
class BoundingBox:
    n: float
    e: float
    s: float
    w: float

    def center(self) -> Point:
        """Find the center of the bounding box.

        Take into account that the bounding box might wrap around the 180°
        meridian.

        >>> BoundingBox(1, 1, -1, -1).center()
        Point(longitude=0.0, latitude=0.0)
        >>> BoundingBox(1, -179, -1, 177).center()
        Point(longitude=179.0, latitude=0.0)
        >>> BoundingBox(1, -177, -1, 179).center()
        Point(longitude=-179.0, latitude=0.0)

        """
        ew = ((self.e + (self.w if self.w < self.e else self.w + 360)) / 2. + 180.) % 360. - 180.
        return Point(ew, (self.n + self.s)/2.)

    def contains(self, p: Point) -> bool:
        """Check whether the bounding box contains some point.

        Take into account that the bounding box might wrap around the 180°
        meridian.

        >>> b = BoundingBox(1, 1, -1, -1)
        >>> b.contains(b.center())
        True
        >>> b = BoundingBox(1, -179, -1, 177)
        >>> b.contains(b.center())
        True
        >>> b = BoundingBox(1, -177, -1, 179)
        >>> b.contains(b.center())
        True
        """

        if self.w < self.e:
            return self.w < p.longitude < self.e and self.s < p.latitude < self.n
        else:
            return (self.e > p.longitude or self.w < p.longitude) and self.s < p.latitude < self.n


land_shp_fname = shpreader.natural_earth(
    resolution='50m', category='physical', name='land')

land_geom = unary_union(
    [record.geometry
     for record in shpreader.Reader(land_shp_fname).records()
     if record.attributes.get('featurecla') != "Null island"])

LAND = prep(land_geom)

AMERICAS = BoundingBox(
    n=74.526716,
    e=-34.535395,
    s=-56.028198,
    w=-168.571541,
)

def coordinates_to_index(points, resolution=2 * 60):
    """Convert long,lat coordinate pairs into indices in a TIF

    Convert a [..., 2] ndarray, or a pair of coordinates, into the matching
    grid indices of a Mercator projection pixel map with a given resolution in
    pixels per degree.

    Paramaters
    ==========
    points: ndarray-like, shape=(..., 2)
        An array of longitude. latitude pairs to be converted to grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray(int), shape=(..., 2)
        An integer array of grid indices

    """
    points = numpy.asarray(points)
    return numpy.stack(
        (numpy.round((-points[..., 1] + 90) * resolution).astype(int),
         numpy.round((points[..., 0] + 180) * resolution).astype(int)),
        -1)


def index_to_coordinates(indices, resolution=2 * 60):
    """Convert grid indices into long,lat coordinate pairs

    Convert a (ndarray of) grid indices of a Mercator projection pixel map with
    a given resolution in pixels per degree into geocoordinate pairs (longitude,
    latitude).

    Paramaters
    ==========
    indices: ndarray(int), shape=(..., 2)
        An integer array of grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray, shape=(..., 2)
        An array of longitude. latitude pairs

    """
    indices = numpy.asarray(indices)
    return numpy.stack(
        (indices[..., 1] / resolution - 180,
         90 - indices[..., 0] / resolution),
        -1)


def geographical_distance(index1: Index, index2: Index) -> float:
    """Calculate the geodesic distance between two hex centers, in meters.

    """
    lat1, lon1 = h3.h3_to_geo(index1)
    lat2, lon2 = h3.h3_to_geo(index2)
    return numpy.asarray(GEODESIC.inverse(
        (lon1, lat1),
        (lon2, lat2)))[0, 0]


class Ecoregions:
    @staticmethod
    def ecoregions():
        """
        >>> eco = Ecoregions.ecoregions()
        >>> eco.numRecords
        847
        >>> eco.fields
        [('DeletionFlag', 'C', 1, 0), ['OBJECTID', 'N', 32, 10], ['ECO_NAME', 'C', 150, 0], ['BIOME_NUM', 'N', 32, 10], ['BIOME_NAME', 'C', 254, 0], ['REALM', 'C', 254, 0], ['ECO_BIOME_', 'C', 254, 0], ['NNH', 'N', 11, 0], ['ECO_ID', 'N', 11, 0], ['SHAPE_LENG', 'N', 32, 10], ['SHAPE_AREA', 'N', 32, 10], ['NNH_NAME', 'C', 64, 0], ['COLOR', 'C', 7, 0], ['COLOR_BIO', 'C', 7, 0], ['COLOR_NNH', 'C', 7, 0], ['LICENSE', 'C', 64, 0]]

        """
        zipshape = zipfile.ZipFile(
            (Path(__file__).parent /
            "../ecoregions/Ecoregions2017.zip").open("rb"))
        shape = shapefile.Reader(
            shp=zipshape.open("Ecoregions2017.shp"),
            shx=zipshape.open("Ecoregions2017.shx"),
            dbf=zipshape.open("Ecoregions2017.dbf"),
            encoding='latin-1'
        )
        return shape

    def __init__(self, mask: t.Optional[sgeom.Polygon] = None):
        shp = self.ecoregions()
        self.boundaries = []
        self.record_getter = shp.record
        for b, boundary in enumerate(shp.shapes()):
            if mask is None:
                self.boundaries.append(prep(sgeom.shape(boundary)))
            else:
                if b != 205 and sgeom.box(*boundary.bbox).intersects(mask):
                    print(self.record_getter(b))
                    s = sgeom.shape(boundary).buffer(0)
                    self.boundaries.append(prep(s & mask))
                else:
                    self.boundaries.append(None)

    def at_point(self, point):
        """Get the ecoregion for the given location

        >>> ECOREGION = Ecoregions()
        >>> r = ECOREGION.at_point(sgeom.Point(-38.282, -8.334))
        >>> r
        84
        >>> ECOREGION.record_getter(r)
        Record #84: [86.0, 'Caatinga', 2.0, 'Tropical & Subtropical Dry Broadleaf Forests', 'Neotropic', 'NO02', 4, 525, 113.696454941, 60.155581182, 'Nature Imperiled', '#81B50A', '#CCCD65', '#EE1E23', 'CC-BY 4.0']
        """
        for b, boundary in enumerate(self.boundaries):
            if boundary and boundary.contains(point):
                return b
        return None

def db():
    engine = sqlalchemy.create_engine("sqlite:///hexbins.sqlite")
    metadata = sqlalchemy.MetaData()
    rectangular_grid_data = sqlalchemy.Table(
        'rect', metadata,
        sqlalchemy.Column('file', sqlalchemy.String, primary_key=True),
        sqlalchemy.Column('x', sqlalchemy.Integer, primary_key=True),
        sqlalchemy.Column('y', sqlalchemy.Integer, primary_key=True),
        sqlalchemy.Column('latitude', sqlalchemy.Float),
        sqlalchemy.Column('longitude', sqlalchemy.Float),
        sqlalchemy.Column('land', sqlalchemy.Boolean),
        sqlalchemy.Column('elevation', sqlalchemy.Integer),
        sqlalchemy.Column('biome', sqlalchemy.Integer),
        sqlalchemy.Column('hexbin', sqlalchemy.Integer),
    )
    hexagonal_grid_data = sqlalchemy.Table(
        'hex', metadata,
        sqlalchemy.Column('hexbin', sqlalchemy.Integer),
    )
    hexagonal_grid_data = sqlalchemy.Table(
        'hex_distance', metadata,
        sqlalchemy.Column('hexbin1', sqlalchemy.Integer),
        sqlalchemy.Column('hexbin1', sqlalchemy.Integer),
        sqlalchemy.Column('distance', sqlalchemy.Float),
    )
    metadata.create_all(engine)



def gmted_tile_from_geocoordinates(lon, lat):
    """
    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010N30E000/
    >>> gmted_tile_from_geocoordinates(15, 40)
    ('GMTED2010N30E000', '30n000e')

    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010S30W090/
    >>> gmted_tile_from_geocoordinates(-75, -20)
    ('GMTED2010S30W090', '30s090w')
    """
    southwest_corner_lon = int(lon // 30) * 30
    southwest_corner_lat = int((lat - 10) // 20) * 20 + 10
    ew = "W" if southwest_corner_lon < 0 else "E"
    ns = "S" if southwest_corner_lat < 0 else "N"
    tile = "GMTED2010{:1s}{:02d}{:1s}{:03d}".format(
        ns, abs(southwest_corner_lat), ew, abs(southwest_corner_lon))
    file = "{:02d}{:1s}{:03d}{:1s}".format(
        abs(southwest_corner_lat), ns.lower(),
        abs(southwest_corner_lon), ew.lower())
    return tile, file

from heapq import heappush as push, heappop as pop
from itertools import count


def navigation_speed(slope):
    """Using the slope in %, calculate the navigation speed in m/s

    This function calculates the off-road navigation speed (for male cadets in
    forested areas with navigational targets every 20 minutes) following
    [@irmischer2018measuring].

    > [T]he fastest off-road navigation speed was 0.78 m/s for males […] with a
    > peak at −2%.

    >>> navigation_speed(-2.)
    0.78
    >>> navigation_speed(-2.) > navigation_speed(-1.5)
    True
    >>> navigation_speed(-2.) > navigation_speed(-2.5)
    True

    """
    ...
    return 0.11 + 0.67 * numpy.exp(-(slope + 2.0) ** 2 / 1800.)

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

eco_cache = {}

def geo_dist(lon0, lat0, lon1, lat1, cache={}):
    try:
        return cache[abs(int(lat0 * 3600)),
                     abs(int(lat1 * 3600)),
                     abs(int((lon1-lon0) * 3600))]
    except KeyError:
        d = GEODESIC.inverse((lon0, lat0), (lon1, lat1))[0, 0]
        cache[abs(int(lat0 * 3600)),
              abs(int(lat1 * 3600)),
              abs(int((lon1-lon0) * 3600))]= d
        return d

def run_on_one_tile(p):
    lat, lon = h3.h3_to_geo(p)
    tile, file = gmted_tile_from_geocoordinates(lon, lat)
    gdal_path = f"/vsizip/../elevation/GMTED2010/{tile:}_150.zip/{file:}_20101117_gmted_med150.tif"
    with rasterio.open(gdal_path) as src:
        ecology = eco_cache[gdal_path] = Ecoregions(sgeom.box(*src.bounds))

        return areas(src, ecology)



def areas(
        raster: rasterio.DatasetReader,
        ecoregions: Ecoregions,
) -> t.Tuple[t.Dict[h3.H3Index, t.Dict[h3.H3Index, float]],
             t.Dict[h3.H3Index, t.Counter[int]]]:
    d = {}
    e = {}
    i, j = raster.shape
    lon, lat = raster.transform * (i/2, j/2)
    starts = {h3.geo_to_h3(lat, lon, RESOLUTION)}
    failed = set()
    data = raster.read(1)
    while starts:
        start = starts.pop()
        print(start)
        try:
            neighbors, metadata = distances_from_focus(
                start,
                data,
                raster,
                ecoregions)
            d[start] = neighbors
            e[start] = metadata
            for n in neighbors:
                if n in d:
                    continue
                starts.add(n)

        except (IndexError, NotImplementedError):
            print("failed")
            failed.add(start)
    return d, e


def distances_from_focus(
        i: h3.H3Index,
        data: numpy.array,
        raster: rasterio.DatasetReader,
        ecoregions: Ecoregions,
) -> t.Tuple[t.Dict[h3.H3Index, float],
             t.Counter[int]]:
        neighbors: t.List[h3.H3Index] = list(h3.k_ring(i, 2))
        origin = neighbors.index(i)
        points = numpy.array([h3.h3_to_geo(n)[::-1] for n in neighbors])
        indices = numpy.round(~raster.transform * points.T).astype(int).T
        source = tuple(indices[origin])
        destinations = {tuple(j) for j in indices}

        # Dijkstra's algorithm, adapted for our purposes from
        # networkx/algorithms/shortest_paths/weighted.html
        dist: t.Dict[t.Tuple[int, int], float] = {}  # dictionary of final distances
        seen = {source: 0}
        c = count()
        fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []  # use heapq with (distance, label) tuples
        push(fringe, (0, next(c), source))


        e = t.Counter()
        def moore_neighbors(x0, y0, biome):
            lon0, lat0 = raster.transform * (x0, y0)
            elevation = data[y0, x0]
            for i, j in [(-1, -1), (-1, 0), (-1, 1),
                         (0, -1), (0, 1),
                         (1, -1), (1, 0), (1, 1)]:
                x1, y1 = x0 + i, y0 + j
                lon1, lat1 = raster.transform * (x1, y1)
                horizontal = geo_dist(lon0, lat0, lon1, lat1)
                ed = data[y1, x1] - elevation
                yield (x1, y1), horizontal / navigation_speed(ed / horizontal * 100) * terrain_coefficients[biome]

        while fringe:
            (d, _, (x0, y0)) = pop(fringe)
            if (x0, y0) in dist:
                continue  # already searched this node.
            dist[(x0, y0)] = d
            if (x0, y0) in destinations:
                destinations.remove((x0, y0))
                if not destinations:
                    break

            lon0, lat0 = raster.transform * (x0, y0)
            j = h3.geo_to_h3(lat0, lon0, RESOLUTION)
            eco = ecoregions.at_point(sgeom.Point(lon0, lat0))
            if j == i:
                e[eco] += 1

            if eco is None:
                print("No biome found at", lon0, lat0)
                continue

            biome = ecoregions.record_getter(eco)[3]


            for u, cost in moore_neighbors(x0, y0, biome):
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
        return {n: dist[x, y] for n, (x, y) in zip(neighbors, indices)}, e

run_on_one_tile('8581a90ffffffff')

def plot_distances(dist):
    r = min(dist, key=dist.get)
    lat0, lon0 = h3.h3_to_geo(r)
    for key, d in dist.items():
        lat1, lon1 = h3.h3_to_geo(key)
        plt.plot((lon0, lon1), (lat0, lat1), linewidth=(200000/(d+1))**2, c='b')

def hexbin_data(start):
        ins = rectangular_grid_data.insert()
        for _, window in raster.block_windows(1):
            r = raster.read(1, window=window)
            trafo = rasterio.windows.transform(window, raster.transform)
            indices = numpy.unravel_index(numpy.arange(len(r.flat)), r.shape)
            lons, lats = trafo * (indices[1], indices[0])
            data = []
            for lon, lat, x, y, value in zip(
                    lons, lats,
                    indices[0] + window.row_off,
                    indices[1] + window.col_off,
                    r.flat):
                index = h3.geo_to_h3(lat, lon, RESOLUTION)
                point = sgeom.Point(lon, lat)
                if value < 10:
                    land = LAND.contains(point)
                else:
                    land = True
                if not land:
                    biome = None
                else:
                    biome = ECOREGION.at_point(point)
                data.append({
                    "file": file,
                    "x": x,
                    "y": y,
                    "latitude": lat,
                    "longitude": lon,
                    "land": land,
                    "biome": biome,
                    "elevation": value,
                    "hexbin": index})
            print(file, x, y)
            engine.execute(ins, data)
