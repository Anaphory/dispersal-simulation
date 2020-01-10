from popcaps import random_popcap
import bisect
import numpy

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
        grow_into = {}
        growth = 0
        for cell in self.cells:
            if cell not in grow_into:
                grow_into[cell] = cell.popcap - cell.population
            region_popcap = cell.popcap
            region_population = cell.population
            for n in cell.neighbors():
                if n not in grow_into:
                    grow_into[n] = n.popcap - n.population
                region_popcap += n.popcap
                region_population += n.population
            if region_popcap == 0:
                continue
            growth += self.growth_factor * cell.population * (1 - region_population / region_popcap)
        self.distribute(growth, list(grow_into))

    def distribute(self, growth, candidates):
        print(growth)
        if not set(candidates) - self.cells:
            raise StopIteration
        if growth < 1:
            raise StopIteration

        for _ in range(int(growth)):
            capacity = [max(1 - cell.population / cell.popcap , 0)
                        for cell in candidates]
            # The paper states distribution is according to a function of N/K,
            # without any specification of the process.
            c = numpy.cumsum(capacity)
            self.popcap -= 1
            if self.popcap <= 0:
                raise StopIteration
            i = bisect.bisect(c, numpy.random.random()*c[0])
            candidates[i].population += 1
            candidates[i].language = self
            self.cells.add(candidates[i])


class FastSpreadLanguage(Language):
    def distribute(self, growth, candidates):
        print(growth)
        if not set(candidates) - self.cells:
            raise StopIteration
        if growth < 1:
            raise StopIteration

        capacities = [max(1 - cell.population / cell.popcap , 0)
                    for cell in candidates]
        # The paper states distribution is according to a function of N/K,
        # without any specification. That is very underspecified and unlikely, though.
        c = sum(capacities)
        for candidate, capacity in zip(candidates, capacities):
            candidate.population += capacity * (growth / c)
            self.popcap -= capacity * (growth / c)
            candidate.language = self
            self.cells.add(candidate)
        if self.popcap <= 0:
            raise StopIteration


class DifferenceSpreadLanguage(Language):
    def distribute(self, growth, candidates):
        print(growth)
        if not set(candidates) - self.cells:
            raise StopIteration
        if growth < 1:
            raise StopIteration

        capacities = [max(cell.popcap - cell.population, 0)
                    for cell in candidates]
        # The paper states distribution is according to a function of N/K,
        # without any specification. That is very underspecified and unlikely, though.
        c = sum(capacities)
        for candidate, capacity in zip(candidates, capacities):
            candidate.population += capacity * (growth / c)
            self.popcap -= capacity * (growth / c)
            candidate.language = self
            self.cells.add(candidate)
        if self.popcap <= 0:
            raise StopIteration


class DifferenceSpreadAndRoundLanguage(Language):
    """This implementation seems to best capture the spirit of the reference.

    It is far less fancy for handling neighborhoods separately and it ignores
    several ways to obtain rounding errors, compared to the original.

    """
    def distribute(self, growth, candidates):
        print(growth)
        if not set(candidates) - self.cells:
            raise StopIteration("no new cells")
        if growth < 1:
            raise StopIteration("no growth")

        sum_growth_space = sum([cell.popcap - cell.population for cell in candidates])
        cell_indices = list(range(len(candidates)))
        while cell_indices:
            i = cell_indices.pop(numpy.random.randint(len(cell_indices)))
            cell = candidates[i]
            cell_growth = round(growth * (cell.popcap - cell.population) / sum_growth_space)
            if cell_growth > 0 or cell.popcap < 1:
                cell.population += cell_growth
                cell.language = self
                self.cells.add(cell)
                print(cell, cell.population, cell_growth)
                self.popcap -= cell_growth
                if self.popcap <= 0:
                    raise StopIteration("popcap reached")
                growth -= cell_growth
            if growth < 0.5:
                return
        while candidates and growth > 1:
            i = numpy.random.randint(len(candidates))
            cell = candidates[i]
            if cell.popcap > cell.population:
                cell.population += 1
                growth -= 1
                self.popcap -= 1
            else:
                del candidates[i]

