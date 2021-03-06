# cython: language_level=3
"""Helper functions for the dispersal model"""

import os
import random

import numpy
import cython
import tifffile

import dispersal_model.osm as osm
from dispersal_model.hexgrid import AMERICAS

from dispersal_model.types_and_units import (
    kcal,
    Any, Callable, Dict, Iterable, Mapping, S, Sequence, T, Tuple, Type)


class PopulationCapModel:
    def __init__(self, alpha: float = 10 ** -8.07, beta: float = 2.64):
        self.alpha = alpha
        self.beta = beta
        self.resolution = 60 // 5  # 5 arc minutes = 1/12 degree

    def coordinates_to_index(self, points: Tuple[float, float]) -> Tuple[int, int]:
        """Convert long,lat coordinate pairs into indices in a TIF

        Convert a [..., 2] ndarray, or a pair of coordinates, into the matching
        grid indices of a Mercator projection pixel map with a given resolution
        in pixels per degree.

        Paramaters
        ==========
        points: ndarray-like, shape=(..., 2)
            Array of longitude/latitude pairs to be converted to grid indices

        resolution:
            The resolution of the grid in indices per degree

        Returns
        =======
        ndarray(int), shape=(..., 2)
            An integer array of grid indices

        """
        points_array = numpy.asarray(points)
        return numpy.stack(
            (numpy.round((-points_array[..., 1] + 90) *
                         self.resolution).astype(int),
             numpy.round((points_array[..., 0] + 180) *
                         self.resolution).astype(int)),
            -1)

    @property
    def precipitation_tif(self) -> numpy.ndarray:
        try:
            return self._tif
        except AttributeError:
            self._tif: numpy.ndarray = tifffile.imread(os.path.join(
                os.path.dirname(__file__),
                "wc2.1_5m_bio_12.tif"))
            # WorldClim 2.1 Bioclimate data from Fick & Hijmans (2017), at
            # https://worldclim.org/data/worldclim21.html where Bioclimatic
            # variable 12 is Annual Precipitation according to
            # https://worldclim.org/data/bioclim.html. The 5-minute-resolution
            # data should be slightly finer than our hexgrid in general (and in
            # particular away from the equator). It can be downloaded from
            # https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_5m_bio.zip
            self._tif.clip(0, out=self._tif)
            return self._tif

    def population_capacity(self, point: Tuple[float, float]) -> kcal:
        """Calculate the pop cap of a cell given its precipitation

        Return the carrying capacity K of a hexagonal cell of area AREA with
        the given mean yearly precipitation measured in mm

        In Gavin et al. 2017, the function used is

        K = α * P ** β

        with eg. α = 10 ** -8.07, β = 2.64 [Note: Their supplementary material
        lists α=-8.07, but that is a nonsensical number. The actual plot line
        does not exactly correspond to α=10^-8.07, but more closely to
        α=10^-7.96, but that suggests that this is at least close to the
        described behaviour.]

        Parameters
        ==========
        precipitation: float
            The cell's mean yearly precipitation P, in mm

        Returns
        =======
        float
            The cell's carrying capacity, in individuals/km^2
        """
        return self.alpha * self.precipitation(point) ** self.beta

    def precipitation(self, point: Tuple[float, float]) -> int:
        index = tuple(self.coordinates_to_index(point))
        return self.precipitation_tif[index]


p = PopulationCapModel(
    alpha=4 * 10 ** -8.07,
    # First inspection suggest an under-estimation of the population of Alaska
    # by a factor of ~4, so we add that factor here for now. The goal is to
    # provide a better population model in the long run, so this construction
    # will hopefully be completely replaced anyway.
    beta=2.64)


@cython.ccall
def density(land: bool, longitude: float, latitude: float) -> float:
    """Estimate the population density in person/km²."""
    if longitude < AMERICAS.w or longitude > AMERICAS.e:
        return 0
    if latitude < AMERICAS.s or latitude > AMERICAS.n:
        return 0
    if not land:
        return 0
    else:
        return p.population_capacity(
            (longitude, latitude))


def shuffle(seq: Sequence) -> Iterable:
    """Return a shuffled version of `seq`.

    seq is copied, not modified, in the process.

    >>> w = [1, 3, 2, 4]
    >>> sorted(set(shuffle(w)))
    [1, 2, 3, 4]
    >>> w
    [1, 3, 2, 4]

    """
    s = list(seq)[:]
    while s:
        i = random.randrange(len(s))
        yield s[i]
        del s[i]


D: Type[dict]
if cython.compiled:
    D = dict
else:
    D = Dict[S, T]


class OnDemandDict(D):
    """A dictionary whose missing keys are filled on-demand from a Callable.

    >>> d = OnDemandDict(lambda key: key)
    >>> d[1]
    1
    >>> d[2]
    2
    >>> d
    {1: 1, 2: 2}

    """
    def __init__(self, populating_function: Callable[[S], T], **kwargs):
        super().__init__(**kwargs)
        self.function = populating_function

    def __missing__(self, key: S) -> T:
        self[key] = self.function(key)
        return self[key]


def serialize(o: object) -> object:
    """JSON helper function to serialize nummpy integers"""
    if isinstance(o, numpy.integer):
        return int(o)
    raise TypeError


def get_data(longitude: float, latitude: float) -> Mapping[str, Any]:
    "Get climate data for this point"
    return {
        "longitude": longitude,
        "latitude": latitude,
        "land": osm.is_land((longitude, latitude))}


def random_choice(x: Sequence[T]) -> T:
    """Return a random element from the sequence."""
    return x[random.randrange(len(x))]


def in_random_order_ignoring_location(
        keyval: Mapping[Any, Sequence[T]]) -> Iterable[T]:
    """Return the members of the mappings's values in random order.

    >>> set(in_random_order_ignoring_location(
    ... {1: [2, 3], 4: [5], 6: [], 7: [1]}))
    {1, 2, 3, 5}

    """
    return shuffle(sum((f for f in keyval.values()), []))

