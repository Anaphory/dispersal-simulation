"""Benchmark different ways of calculating culture similarity.

Provide the fastest one under a generic name for imports.

"""

import cython


@cython.inline
@cython.ccall
def cultural_distance(c1: cython.ulong, c2: cython.ulong) -> cython.ulong:
    """Cultural distance is the Hamming distance of the culture vectors.

    Because cultures are stored in binary, the Hamming distance is all the bits
    in c1 XOR c2. This may look a bit intransparent to the non-programmer, but
    it works:

    >>> cultural_distance(0b001000, 0b000101)
    3
    >>> cultural_distance(0b1110111, 0b1110111)
    0
    >>> cultural_distance(0b011011011, 0b010011011)
    1

    Obviously, this also works with the more compact decimal representation of
    binary vectors – i.e. numbers.

    >>> 8 == 0b1000
    True
    >>> 3 == 0b0011
    True
    >>> cultural_distance(8, 3) == 3 == cultural_distance(0b0100, 0b0011)
    True

    """
    return bin(c1 ^ c2).count('1')


@cython.inline
@cython.ccall
def similar_culture_1(c1: cython.ulong, c2: cython.ulong) -> int:  # bool
    return cultural_distance(c1, c2) < 6


@cython.inline
@cython.ccall
def similar_culture_2(c1: cython.ulong, c2: cython.ulong) -> int:  # bool
    return bin(c1 ^ c2).count('1') < 6


@cython.inline
@cython.ccall
def similar_culture_3(c1: cython.ulong, c2: cython.ulong) -> int:  # bool
    f = c1 ^ c2
    return (
        (f & 1) +
        (f >> 1 & 1) +
        (f >> 2 & 1) +
        (f >> 3 & 1) +
        (f >> 4 & 1) +
        (f >> 5 & 1) +
        (f >> 6 & 1) +
        (f >> 7 & 1) +
        (f >> 8 & 1) +
        (f >> 9 & 1) +
        (f >> 10 & 1) +
        (f >> 11 & 1) +
        (f >> 12 & 1) +
        (f >> 13 & 1) +
        (f >> 14 & 1) +
        (f >> 15 & 1) +
        (f >> 16 & 1) +
        (f >> 17 & 1) +
        (f >> 18 & 1) +
        (f >> 19 & 1))


@cython.inline
@cython.ccall
def similar_culture_4(c1: cython.ulong, c2: cython.ulong) -> int:  # bool
    c = 0
    f = c1 ^ c2
    while f and c < 6:
        c += f & 1
        f >>= 1
    return c < 6


similar_culture = similar_culture_1

if __name__ == "__main__":
    # Run the benchmark
    import timeit
    for i in range(1, 5):
        print(
            f"similar_culture_{i:}",
            timeit.timeit(
                f"""for c in range(2<<20):
  similar_culture_{i:}(0b0101010101010101010101010101010101010101, c)""",
                number=1,
                globals=globals()))
