from popcaps import random_popcap
import bisect
import numpy
from collections import Counter

class Language:
    growth_factor = 1.1

    def __init__(self, color, cell):
        g = cell
        g.population = 10
        g.language = self
        self.cells = {g}
        self.id = color
        self.popcap = random_popcap()

    def grow(self):
        grow_into = Counter()
        for cell in self.cells:
            region_growth = set()
            region_growth.add(cell)
            region_popcap = cell.popcap
            region_population = cell.population

            for n in cell.neighbors():
                region_growth.add(n)
                region_popcap += n.popcap
                region_population += n.population

            if region_popcap < 1:
                continue

            growth = self.growth_factor * cell.population * (1 - region_population / region_popcap)
            if growth < 0:
                continue

            grow_into += self.distribute(growth, region_growth)

        for cell, growth in grow_into.items():
            if cell.popcap < 1 or growth > 1:
                cell.population += growth
                cell.language = self
                self.cells.add(cell)

    def distribute(self, growth, candidates):
        """This implementation seems to best capture the spirit of the reference.

        It is far less fancy for handling neighborhoods separately and it ignores
        several ways to obtain rounding errors, compared to the original.

        """
        print(growth)

        sum_growth_space = sum([cell.popcap - cell.population for cell in candidates])
        
        cell_growth = Counter()

        cell_indices = [0] + sorted(range(1, len(candidates)), key=lambda _: numpy.random.random())
        for i in cell_indices:
            cell = list(candidates)[i]
            cell_growth[cell] = int(
                min(growth, cell.popcap - cell.population))
            print(cell, cell.population, cell_growth[cell])
            self.popcap -= cell_growth[cell]
            if self.popcap <= 0:
                raise StopIteration("popcap reached")
            growth -= cell_growth[cell]

        return cell_growth


class DifferenceSpreadAndRoundLanguage(Language):
    pass
