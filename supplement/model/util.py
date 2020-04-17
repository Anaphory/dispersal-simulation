"""Helper functions for the dispersal model"""

import numpy

import osm
from hexgrid import AMERICAS

from types_and_units import (
    Sequence,
    Dict,
    Iterable,
    S,
    T,
    Callable,
    Mapping,
    Iterator,
    Any,
)


from gavin2017processbased import PopulationCapModel, Point
p = PopulationCapModel()
p.alpha = 4 * 10 ** -8.07
p.beta = 2.64


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
            Point(longitude=longitude, latitude=latitude))


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


def in_random_order_ignoring_location(keyval: Mapping[Any, Sequence[T]]) -> Iterable[T]:
    """Return the members of the mappings's values in random order.

    >>> set(in_random_order_ignoring_location({1: [2, 3], 4: [5], 6: [], 7: [1]}))
    {1, 2, 3, 5}

    """
    return shuffle(sum((f for f in keyval.values()), []))
