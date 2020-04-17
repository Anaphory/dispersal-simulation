"""Helper functions for the dispersal model"""

import numpy

import osm
import tifffile
from hexgrid import AMERICAS

from types_and_units import (
    Sequence,
    Dict,
    Iterable,
    S,
    T,
    Callable,
    Mapping,
    Any,
)


class PopulationCapModel:
    def __init__(self, alpha: float = 10 ** -8.07, beta: float = 2.64):
        self.alpha = alpha
        self.beta = beta
        self.resolution = 2 * 60

    def coordinates_to_index(self, points):
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
        points = numpy.asarray(points)
        return numpy.stack(
            (numpy.round((-points[..., 1] + 90) *
                         self.resolution).astype(int),
             numpy.round((points[..., 0] + 180) *
                         self.resolution).astype(int)),
            -1)

    @property
    def precipitation_tif(self):
        try:
            return self._tif
        except AttributeError:
            self._tif = tifffile.imread("/home/gereon/Public/settlement-of-americas/supplement/worldclim/wc2.0_bio_30s_12.tif").clip(0)
            return self._tif

    def population_capacity(self, point):
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

    def precipitation(self, point):
        index = tuple(self.coordinates_to_index(point))
        return self.precipitation_tif[index]


p = PopulationCapModel(
    alpha=4 * 10 ** -8.07,
    beta=2.64)


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
            [longitude, latitude])


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
        i = numpy.random.randint(len(s))
        yield s[i]
        del s[i]


class OnDemandDict(Dict[S, T]):
    """A dictionary whose missing keys are filled on-demand from a Callable.

    >>> d = OnDemandDict(lambda key: key)
    >>> d[1]
    1
    >>> d[2]
    2
    >>> d
    {1: 1, 2: 2}

    """
    def __init__(self, populating_function: Callable[[S], T], *args, **kwargs):
        super().__init__(*args, **kwargs)
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
    return x[numpy.random.randint(len(x))]


def in_random_order_ignoring_location(
        keyval: Mapping[Any, Sequence[T]]) -> Iterable[T]:
    """Return the members of the mappings's values in random order.

    >>> set(in_random_order_ignoring_location(
    ... {1: [2, 3], 4: [5], 6: [], 7: [1]}))
    {1, 2, 3, 5}

    """
    return shuffle(sum((f for f in keyval.values()), []))
