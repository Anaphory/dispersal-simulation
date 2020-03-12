import sys
import attr
import math
import numpy
import pickle
import argparse
import itertools
from abc import abstractmethod
import cartopy.crs as ccrs
import gavin2017processbased
import matplotlib
from matplotlib import pyplot as plt
from gavin2017processbased import is_land

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

coordinates: numpy.ndarray = gavin2017processbased.hexagonal_earth_grid(
    gavin2017processbased.americas,
    gavin2017processbased.area)

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

def hexagon_coords(i: Index, neighbors: Sequence[Index]):
    n = iter(neighbors)
    old = first = next(n)
    l = []
    for j in n:
        l.append(list((coordinates[j] + coordinates[old] + coordinates[i]) / 3))
        old = j
    l.append(list((coordinates[first] + coordinates[old] + coordinates[i]) / 3))
    return l

def plot(family_locations, hexes, maximum=None) -> None:
    plt.gcf().set_size_inches(30, 30)
    polygons = []
    values = []
    for cell_coords, cell_value in hexes:
        polygons.append(numpy.array(cell_coords))
        values.append(cell_value)
    cmap = plt.get_cmap("viridis")
    vmax = maximum or max(values)
    values = [cmap(255 * v / vmax) for v in values]
    collection = matplotlib.collections.PolyCollection(
        polygons, facecolors=values)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_collection(collection)
    ax.coastlines("50m")
    ax.set_extent(gavin2017processbased.americas)

    coords = []
    for family, location in family_locations.items():
        m, i, j = location
        coords.append(coordinates[m][i, j] + numpy.random.normal(0, 0.1, 2))
    plt.scatter(*zip(*coords), alpha=0.3, s=3, edgecolors=None)

def plot_(s):
    plot({family.descendence: family.location for family in s.families},
         [(hexagon_coords(index), patch.resources)
          for index, patch in s.grid.patches.items()])

    plt.show()

def plot_last_line(filename):
    with open(filename, "r") as f:
        for line in f:
            pass
    properties = eval(line.strip("\x00\n"))
    # DAFUQ? Why even are there 00 bytes??
    families = properties["Locations"]
    hexes = properties["Resources"]
    plot(families, hexes)
    plt.show()

def plot_series(filename, template="dispersal-{:04d}.png"):
    with open(filename, "r") as f:
        for l, line in enumerate(f):
            properties = eval(line.strip("\x00\n"))
            families = properties["Locations"]
            hexes = properties["Resources"]
            plot(families, hexes)
            plt.savefig(template.format(l))
            plt.close()


def geographical_distance(index1: Index, index2: Index) -> meters:
    m1, i1, j1 = index1
    m2, i2, j2 = index2
    return numpy.asarray(GEODESIC.inverse(
        coordinates[m1][i1, j1], coordinates[m2][i2, j2])[:, 0])[0]

def neighbors_within_distance(coordinates: Tuple[numpy.ndarray, numpy.ndarray],
                              mij: Index, d: meters) -> Iterable[Index]:
    """A generator of all cells within d meters of mij.

    Yield mij back as first item of the generator.
    """
    aggregate: List[Index] = []
    m, i, j = mij
    yield m, i, j
    aggregate.append(mij)
    c = coordinates[m][i, j]
    for r_m, block in enumerate(coordinates):
        b = block.reshape((-1, 2))
        r_ij = numpy.arange(b.shape[0])[numpy.asarray(GEODESIC.inverse(c, b)[:, 0]) <= d]
        for r_i, r_j in numpy.array(numpy.unravel_index(r_ij, coordinates[r_m].shape[:-1])).T:
            if r_m == m and r_i == i and r_j == j:
                continue
            yield r_m, r_i, r_j
            aggregate.append((r_m, r_i, r_j))
