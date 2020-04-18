"""Typed wrappers around Uber's h3 module for the dispersal model."""

import numpy
import cython
from h3 import h3

import cartopy.geodesic as geodesic

from types_and_units import meters, Iterator, Tuple, NamedTuple

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
