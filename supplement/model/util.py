import sys
import attr
import json
import math
import numpy
import pickle
import argparse
import itertools
import matplotlib
import cartopy.crs as ccrs
from abc import abstractmethod
from matplotlib import pyplot as plt

import osm
import gavin2017processbased
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

def hexagon_coords(i: Index, neighbors: Sequence[Index], cache: Dict[Index, List]={}) -> List:
    try:
        return cache[i]
    except KeyError:
        pass
    n = iter(neighbors)
    old = first = next(n)
    l = []
    for j in n:
        l.append(list((coordinates[j] + coordinates[old] + coordinates[i]) / 3))
        old = j
    l.append(list((coordinates[first] + coordinates[old] + coordinates[i]) / 3))
    cache[i] = l
    return l

def plot(family_locations, hexes, maximum=365000000.0) -> None:
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
    sizes = []
    for family, location in family_locations:
        sizes.append(family)
        coords.append(coordinates[tuple(location)] + numpy.random.normal(0, 0.04, 2))
    plt.scatter(*zip(*coords), alpha=0.2, s=sizes, c="gray", edgecolors=None)

def plot_(s):
    plot([(family.effective_size, family.location) for family in s.families],
         [(hexagon_coords(index), patch.resources)
          for index, patch in s.grid.patches.items()])

    plt.show()

def plot_last_line(filename):
    with open(filename, "r") as f:
        for line in f:
            pass
    properties = json.loads(line)
    # DAFUQ? Why even are there 00 bytes??
    families = properties["Families"].values()
    hexes = properties["Resources"]
    plot(families, hexes)
    plt.show()

def plot_series(filename, template="dispersal-{:07d}.png", limit=None):
    with open(filename, "r") as f:
        for l, line in enumerate(f):
            if limit is not None and l not in limit:
                continue
            properties = json.loads(line)
            families = properties["Families"].values()
            hexes = properties["Resources"]
            plot(families, hexes)
            plt.savefig(template.format(l))
            plt.close()


def plot_alaska_population(filename):
    pop = []
    cache = {}
    def is_in(location):
        try:
            return cache[location]
        except KeyError:
            cache[location] = osm.contains(osm.shapes["Alaska"], coordinates[location])
            return cache[location]
    with open(filename, "r") as f:
        for l, line in enumerate(f):
            pop.append(0)
            try:
                data = json.loads(line)
            except json.JSONDecodeError:
                break
            for effective_size, location in data["Families"].values():
                if is_in(tuple(location)):
                    pop[l] += effective_size
    plt.plot(pop)
    plt.show()

def geographical_distance(index1: Index, index2: Index) -> meters:
    m1, i1, j1 = index1
    m2, i2, j2 = index2
    return numpy.asarray(GEODESIC.inverse(
        coordinates[m1][i1, j1], coordinates[m2][i2, j2])[:, 0])[0]
