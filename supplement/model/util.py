import sys
import attr
import json
import math
import numpy
import pickle
import argparse
import itertools
import matplotlib
import igraph
import cartopy.crs as ccrs
from abc import abstractmethod
from matplotlib import pyplot as plt

import osm
import gavin2017processbased
from gavin2017processbased import is_land

def cultural_distance(c1, c2) -> float:
    """ Cultural distance is the Hamming distance of the culture vectors """
    return sum(e1!=e2 for e1, e2 in zip(c1, c2))

# Import some modules for the Entities, in particular types used to describe
# the Entities
from typing import (
    Any,
    Callable,
    DefaultDict,
    Dict,
    IO,
    Iterable,
    Iterator,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
)
from argparse import Namespace

Index = Tuple
meters = float
kcal = float
halfyears = int


def closest_grid_point(longitude, latitude, coordinates):
    return numpy.unravel_index(numpy.argmin(abs(coordinates - (longitude, latitude)).sum(-1)), coordinates.shape[:-1])

ParametersT = TypeVar("ParametersT", bound="Parameters")
class Parameters (
        NamedTuple("Parameters", [
            ("n_steps", halfyears),
        ])):
    @staticmethod
    def from_ns(ns: Namespace) -> 'Parameters':
        return Parameters(**vars(ns))

def shuffle(seq: Sequence) -> Iterable:
    s = list(seq)[:]
    while s:
        i = numpy.random.randint(len(s))
        yield s[i]
        del s[i]

import cartopy.geodesic as geodesic
GEODESIC = geodesic.Geodesic()

S = TypeVar("S")
T = TypeVar("T")
class OnDemandDict(Dict[S, T]):
    def __init__(self, populating_function: Callable[[S], T], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.function = populating_function
    def __missing__(self, key: S) -> T:
        self[key] = self.function(key)
        return self[key]

def serialize(o):
    if isinstance(o, numpy.integer):
        return int(o)
    raise TypeError

def population_within(area):
    """"""
    ...

def get_data(longitude, latitude):
    "Get climate data for this point"
    return {
        "longitude": longitude,
        "latitude": latitude,
        "land": is_land((longitude, latitude))}

def random_choice(x: Sequence[T]) -> T:
    return x[numpy.random.randint(len(x))]

def cached_generator(function: Callable) -> Callable:
    cache: Mapping[Tuple, Any] = {}
    def func(*args):
        try:
            return cache[tuple(args)]
        except KeyError:
            cache[tuple(args)] = list(function(*args))
            return cache[tuple(args)]
    func.__name__ = function.__name__
    return func

def nearby_cell(latitude, longitude, coordinates):
    return numpy.unravel_index(
        numpy.argmin(
            ((coordinates - numpy.array((longitude, latitude))) ** 2).sum(-1)),
        coordinates.shape[:-1])

def hexagon_coords(i: Index, coordinates:numpy.ndarray, neighbors: Sequence[Index], cache: Dict[Index, List]={}) -> List:
    try:
        return cache[i]
    except KeyError:
        pass
    n = iter(neighbors)
    old = first = next(n)
    sh = coordinates.shape[:-1]
    while (numpy.zeros_like(sh) > old).any() | (old > numpy.array(sh) - 1).any():
        old = first = next(n)
    l = []
    for j in n:
        if (numpy.zeros_like(sh) > j).any() | (j > numpy.array(sh) - 1).any():
            continue
        l.append(list((coordinates[j] + coordinates[old] + coordinates[i]) / 3))
        old = j
    l.append(list((coordinates[first] + coordinates[old] + coordinates[i]) / 3))
    cache[i] = l
    return l

def neighbors(mij: Index):
    m, i, j = mij
    if m==0:
        return [(0, i, j+1), (1, i, j), (1, i, j-1),
                (0, i, j-1), (1, i-1, j-1), (1, i-1, j)]
    else:
        return [(1, i, j+1), (0, i, j+1), (0, i, j),
                (1, i, j-1), (0, i+1, j), (0, i+1, j+1)]

def geographical_distance(index1: Index, index2: Index, coordinates: numpy.ndarray) -> meters:
    m1, i1, j1 = index1
    m2, i2, j2 = index2
    return numpy.asarray(GEODESIC.inverse(
        coordinates[m1][i1, j1], coordinates[m2][i2, j2])[:, 0])[0]

