class GridCell():
    alpha = 10 ** -8.07
    beta = 2.64

    grid = hexagonal_earth_grid(
        BoundingBox(
            -1, 1, -1, 1),
        area)
    all_gridcells = {}

    @classmethod
    def gridcell(k, m, i, j):
        if i < 0 or i >= k.grid[m].shape[0]:
            return None
        if j < 0 or j >= k.grid[m].shape[1]:
            return None
        try:
            return k.all_gridcells[m, i, j]
        except KeyError:
            k.all_gridcells[m, i, j] = k(m, i, j)
            return k.all_gridcells[m, i, j]


    def __init__(self, m, i, j, grid=grid):
        self.m = m
        self.ij = i, j
        self.population = 0
        self.popcap = self.population_capacity() * area / 1000000 # area is in m², popcap is per km²
        self.language = None

    def polygon(self):
        try:
            return self._polygon
        except AttributeError:
            neighbors = numpy.array([n.point for n in self.neighbors(True, True)])
            self._polygon = (neighbors + numpy.roll(neighbors, 1, 0) + self.point) / 3
            return self._polygon

    @property
    def point(self):
        return Point(*self.grid[self.m][self.ij])

    def __hash__(self):
        return hash((self.m, self.ij))

    def __repr__(self):
        return "<Cell {:}:{:},{:}{:}>".format(
            self.m, self.ij[0], self.ij[1],
            " with language {:}".format(self.language.id) if self.language else " (empty)")

    def neighbors(self, include_unlivable=False, include_foreign=False):
        i, j = self.ij
        if self.m==0:
            neighbors = [
                self.gridcell(0, i, j+1),
                self.gridcell(1, i, j),
                self.gridcell(1, i, j-1),
                self.gridcell(0, i, j-1),
                self.gridcell(1, i-1, j-1),
                self.gridcell(1, i-1, j),
            ]
        else:
            neighbors = [
                self.gridcell(1, i, j+1),
                self.gridcell(0, i, j+1),
                self.gridcell(0, i, j),
                self.gridcell(1, i, j-1),
                self.gridcell(0, i+1, j),
                self.gridcell(0, i+1, j+1),
            ]
        neighbors = [g or self for g in neighbors]
        return [g for g in neighbors
                if include_foreign or g.language == self.language or g.language is None
                if include_unlivable or g.popcap >= 1]

    @classmethod
    def random_cell(k):
        m = numpy.random.randint(2)
        return (
            m,
            numpy.random.randint(k.grid[m].shape[0]),
            numpy.random.randint(k.grid[m].shape[1]))

