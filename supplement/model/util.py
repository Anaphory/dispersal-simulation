import sys
import attr
import math
import numpy
import pickle
import argparse

# Import some modules for the Entities, in particular types used to describe
# the Entities
from typing import (
    Mapping,
    Tuple,
    Callable,
    Iterable,
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
    s = seq[:]
    while s:
        i = numpy.random.randint(len(s))
        yield s[i]
        del s[i]

import cartopy.geodesic as geodesic
GEODESIC = geodesic.Geodesic()
def neighbors_within(coordinates):
    def neighbors_within_distance(mij: Index, d: meters) -> Sequence[Index]:
        m, i, j = mij
        c = coordinates[m][i, j]
        for r_m, block in enumerate(coordinates):
            b = block.reshape((-1, 2))
            r_ij = numpy.arange(b.shape[0])[numpy.asarray(GEODESIC.inverse(c, b)[:, 0]) < d]
            for r_i, r_j in numpy.array(numpy.unravel_index(r_ij, coordinates[r_m].shape[:-1])).T:
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
    return {}


