"""drennan2004comparing

In their article, Drennan&Peterson (2004) suggest an ‘A-coefficient’ to measure
the convexity or primateness of a population distribution into subgroups.

Drennan, Robert D & Peterson, Christian E. 2004. Comparing archaeological
settlement systems with rank-size graphs: a measure of shape and statistical
confidence. Journal of Archaeological Science 31(5). 533–549.
(doi:10.1016/j.jas.2003.10.002)

"""
import numpy

def a_coefficient(sequence):
    """Calculate the A coefficient of a distribution of group sizes.

    Examples
    ========
    The A coefficient for a uniform population is 1.

    >>> a_coefficient([3, 3, 3, 3, 3])
    1.0

    For a convex population, the A coefficient is positive.

    >>> a_coefficient([3, 10, 12, 7])
    0.4938544382831819

    In a primate population, with a core group and a bunch of individuals
    outside that group, the A coefficient is negative.

    >>> a_coefficient([1, 1, 9, 1, 1, 1, 1, 1, 1])
    -0.6845351232142713

    For extreme discrepancies, it can even be below -1.

    >>> a_coefficient([1, 1, 999, 1, 1, 1, 1, 1, 1])
    -4.295166972038271

    For a Zipf distribution, where population~1/rank, the A coefficent is 0.
    >>> a_coefficient([1, 1/2, 1/3, 1/4, 1/5])
    0.0

    Parameters
    ==========
    sequence: Sequence
        A sequence of positive numbers

    Returns
    =======
    float
        The A coefficient, a number in (-∞, 1]

    """
    log_sizes = numpy.log(sorted(sequence, reverse=True))
    if not len(log_sizes):
        return numpy.nan
    log_rank = numpy.log(numpy.arange(1, len(log_sizes) + 1))
    area = 0
    z0 = log_sizes[0]
    pr, ps, pz = 0, z0, z0
    area = 0
    for lr, ls in zip(log_rank, log_sizes):
        lz = z0 - lr
        area += (ps - pz) * (lr - pr) + (ls - lz - ps + pz) * (lr - pr) / 2
        pr, ps, pz = lr, ls, lz
    return 2 * area / lr ** 2
