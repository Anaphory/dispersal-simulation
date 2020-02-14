from popcaps import random_popcap
import bisect
import numpy
from collections import Counter

class Language:
    growth_factor = 1.1

    def __init__(self, color, cell):
        g = cell
        g.population = max(1, min(10, int(g.popcap)))
        g.language = self
        self.cells = {g}
        self.id = color
        self.popcap = random_popcap()
        print(self.popcap)

    def grow(self):
        grow_into = Counter()
        growth = 0
        for cell in self.cells:
            region_growth = set()
            region_growth.add(cell)
            region_popcap = cell.popcap
            region_population = cell.population

            for n in cell.neighbors():
                region_growth.add(n)
                region_popcap += n.popcap
                region_population += n.population

            growth += self.growth_factor * cell.population * (1 - region_population / region_popcap)

            more_grow_into, growth = self.distribute(growth, region_growth)
            grow_into += more_grow_into

        grown = False
        for cell, growth in grow_into.items():
                grown = True
                cell.population += growth
                cell.language = self
                self.cells.add(cell)
                self.popcap -= growth
                if self.popcap <= 0:
                    raise StopIteration("popcap reached")
        if not grown:
            raise StopIteration("No more growth")

    def distribute(self, growth, candidates):
        """This implementation seems to best capture the spirit of the reference.

        It is far less fancy for handling neighborhoods separately and it ignores
        several ways to obtain rounding errors, compared to the original.

        """
        growth_space = [cell.popcap - cell.population for cell in candidates]

        cell_growth = Counter()

        cell_indices = [0] + sorted(range(1, len(candidates)), key=lambda _: numpy.random.random())

        m = dhondt(growth, growth_space)
        # Use a divisor method to distribute the right number of new
        # individuals. Divisor methods have a rounding choice. Here we round
        # down (as per D'Hondt), and as such favour the hexes with the most
        # growth potential. Other options would be to round to the nearest
        # integer, or to round up (which favours small areas). The original
        # paper uses a quota method, not a divisor method, with the additional
        # individuals distributed at random, thus slightly favouring small
        # populations, but in a different manner.'.
        for i in cell_indices:
            cell = list(candidates)[i]
            cg = m[i]
            cell_growth[cell] = cg
            growth -= cg
            # Fill cells while their growth is positive, to prevent spilling
            # into empty cells too much, but do walk through desert cells if
            # they by chance appear before other cells still to be filled.
            # (This will not be the case when the focus cell is already a
            # desert cell, because then growth will start at 0.)
            if growth < 1:
                break


        return cell_growth, growth

def dhondt(seats_to_allocate, votes):
    """Allocate seats according to votes using D'Hondt's methods

    In case of a tie, the seat goes to the votes that come earlier in the list.
    """
    d = [1 for v in votes]
    result = [0 for v in votes]
    while sum(result) < seats_to_allocate:
        next_seat = max(range(len(votes)), key=lambda i: votes[i]/d[i])
        result[next_seat] += 1
        d[next_seat] += 1

    return result

class DifferenceSpreadAndRoundLanguage(Language):
    pass
