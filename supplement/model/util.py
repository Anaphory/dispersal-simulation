import sys
import attr
import math
import numpy
import pickle
import argparse
import itertools
from gavin2017processbased import is_land

# Import some modules for the Entities, in particular types used to describe
# the Entities
from typing import (
    Any,
    Mapping,
    Tuple,
    Callable,
    Iterable,
    Iterator,
    IO,
    Sequence,
    DefaultDict,
    NamedTuple,
    TypeVar,
    Optional,
    List,
)
from argparse import Namespace

Index = Tuple
meters = float
kcal = float
halfyears = int

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
def neighbors_within(coordinates):
    cache = {}
    def neighbors_within_distance(mij: Index, d: meters) -> Sequence[Index]:
        """A generator of all cells within d meters of mij.

        Yield mij back as first item of the generator.
        """
        m, i, j = mij
        yield m, i, j
        c = coordinates[m][i, j]
        for r_m, block in enumerate(coordinates):
            b = block.reshape((-1, 2))
            r_ij = numpy.arange(b.shape[0])[numpy.asarray(GEODESIC.inverse(c, b)[:, 0]) <= d]
            for r_i, r_j in numpy.array(numpy.unravel_index(r_ij, coordinates[r_m].shape[:-1])).T:
                if r_m == m and r_i == i and r_j == j:
                    continue
                yield r_m, r_i, r_j
    return neighbors_within_distance


class OnDemandDict(dict):
    def __init__(self, populating_function, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.function = populating_function
    def __missing__(self, key):
        self[key] = self.function(key)
        return self[key]


def get_data(longitude, latitude):
    "Get climate data for this point"
    return {
        "longitude": longitude,
        "latitude": latitude,
        "land": is_land((longitude, latitude))}

T = TypeVar("T")
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
