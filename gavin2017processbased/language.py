from popcaps import random_popcap

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
            growth += self.growth_factor * cell.population * (1 - region_population / region_popcap)
        if not set(grow_into) - self.cells:
            raise StopIteration
        self.popcap -= growth
        if self.popcap < 0:
            raise StopIteration
        distribution = growth / sum(grow_into.values())
        if distribution > 1:
            print("Population will grow by more than its population capacity allows: {:}".format(distribution))
        for cell, proportion in grow_into.items():
            if cell.popcap < 1:
                continue
            cell.language = self
            cell.population += proportion * distribution
        self.cells = self.cells.union(grow_into)
        print(growth)

