import numpy
from typing import (
    Mapping,
    Tuple,
    Callable,
    Iterable,
    Sequence,
    DefaultDict
)

def permit_deletion(seq: Sequence) -> Iterable:
    old = None
    for i in seq:
        yield old
        old = i
    yield old

def shuffle(seq: Sequence) -> Iterable:
    s = seq[:]
    while s:
        i = numpy.random.randint(len(s))
        yield s[i]
        del s[i]
