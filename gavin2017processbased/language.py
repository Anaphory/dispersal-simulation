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
            # without any specification of the process. But it also implicitly
            # states that no randomness is involved in the distribution of new
            # individuals.
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
        # without any specification. This gives beautiful broad regular
        # hexagons instead of their approximate circles, though.
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
        # without any specification. Here I use something based on K-N instead.
        c = sum(capacities)
        for candidate, capacity in zip(candidates, capacities):
            candidate.population += capacity * (growth / c)
            self.popcap -= capacity * (growth / c)
            candidate.language = self
            self.cells.add(candidate)
        if self.popcap <= 0:
            raise StopIteration
