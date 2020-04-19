"""Typed wrappers around Uber's h3 module for the dispersal model."""

import numpy
import cython
from h3 import h3

import cartopy.geodesic as geodesic

from dispersal_model.types_and_units import meters, Iterator, Tuple, NamedTuple

RESOLUTION = 5

AREA = h3.hex_area(RESOLUTION, "km^2")

# The resolution of our hexes should be about 450 km² following Gavin (2017)
assert h3.hex_area(RESOLUTION - 1, "km^2") > 450 > AREA

GEODESIC = geodesic.Geodesic()

Index = str


class BoundingBox(NamedTuple):
    """A bounding box with geo-coordinates west, east, south, and north"""
    w: cython.float
    e: float
    s: float
    n: float


AMERICAS = BoundingBox(
    e=-34.535395,
    s=-56.028198,
    w=-168.571541,
    n=74.526716
)


def neighbors_within_distance(mij: Index, distance: meters) -> Iterator[Index]:
    r = int(distance / h3.edge_length(RESOLUTION, 'm') * 1.5 + 1)
    # FIXME: r here means 'rough', not 'radius' :-P
    for ring in h3.k_ring_distances(mij, r):
        for i in ring:
            yield i


def geo_coordinates(index: Index) -> Tuple[float, float]:
    """Return logitude, latitude"""
    lat, lon = h3.h3_to_geo(index)
    return lon, lat


def closest_grid_point(longitude: float, latitude: float) -> Index:
    return h3.geo_to_h3(latitude, longitude, RESOLUTION)


def geographical_distance(index1: Index, index2: Index) -> meters:
    lat1, lon1 = h3.h3_to_geo(index1)
    lat2, lon2 = h3.h3_to_geo(index2)
    return numpy.asarray(GEODESIC.inverse(
        (lon1, lat1),
        (lon2, lat2)))[0, 0]


# Instead of going through the h3 module, we can copy the relevant code from
# there directly and reduce type conversions.
Index = h3.c_int
ring_size = int(52413 / h3.edge_length(RESOLUTION, 'm') * 1.5 + 1)
array_len = h3.libh3.maxKringSize(ring_size)
KringArray = h3.H3Index * array_len
DistanceArray = h3.c_int * array_len
krings = KringArray()
distances = DistanceArray()


def neighbors_within_distance(address: Index, _=None) -> Iterator[Index]:
    """Get K-Rings for a given hexagon properly split by ring"""
    # Initializes to zeroes by default, don't need to force
    h3.libh3.kRingDistances(address, ring_size, krings, distances)
    for i in range(0, array_len):
        if krings[i] != 0:
            yield krings[i]


def geo_coordinates(address: Index) -> Tuple[float, float]:
    geo_coord = h3.GeoCoord()
    h3.libh3.h3ToGeo(address, h3.byref(geo_coord))
    return (
        h3.mercator_lng(h3.rads_to_degs(geo_coord.lng)),
        h3.mercator_lat(h3.rads_to_degs(geo_coord.lat)))


def closest_grid_point(longitude: float, latitude: float) -> Index:
    geo_coord = h3.GeoCoord(
        h3.degs_to_rads(latitude), h3.degs_to_rads(longitude))
    return h3.libh3.geoToH3(h3.byref(geo_coord), RESOLUTION)
