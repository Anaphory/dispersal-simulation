import json
import numpy
import collections
import cartopy.geodesic as geodesic
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.collections
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
import tifffile
import matplotlib.transforms as mtransforms
from geo_plot import plot_hex_grid

GEODESIC = geodesic.Geodesic()
SQRT3 = 3 ** 0.5

land_shp_fname = shpreader.natural_earth(
    resolution='50m', category='physical', name='land')

land_geom = unary_union(
    [record.geometry
     for record in shpreader.Reader(land_shp_fname).records()
     if record.attributes.get('featurecla') != "Null island"])

LAND = prep(land_geom)


# Define classes
BoundingBox = collections.namedtuple("boundingbox",
    ["w", "e", "s", "n"])
Point = collections.namedtuple("point",
    ["longitude", "latitude"])


def hexagon(area):
    """Calculate the edge length and midpoint distance of a plane hexagon.

    Given a hexagon of area `area`, return its edge length and the distance
    between closest neighbour mid points in a hexagonal grid.

    Parameters
    ==========
    area: number
        The area of a plane hexagon

    Returns
    =======
    hexagon_side: float
        The edge length of the hexagon
    grid_point_distance: float
        The distance beween the mid points of two adjacent hexagons in a hexagonal tiling
    """
    hexagon_side = (2 * SQRT3 * area) ** 0.5 / 3
    grid_point_distance = SQRT3 * hexagon_side
    assert abs((grid_point_distance * hexagon_side * 1.5) - area) < 1e-2
    return hexagon_side, grid_point_distance


def hexagonal_earth_grid(bbox, area):
    """Generate a hexagonal grid.

    Generate a spherical hexagonal grid inside the bounding box in which each
    hexagon has roughly the given average area.

    In the return value, mid_grid[k, l] has neighbours mid_grid[k-1, l],
    mid_grid[k+1, l], grid[k, l], grid[k+1, l], grid[k, l+1], grid[k+1, l+1]

    Parameters
    ==========
    area: float
        The desired hexagon surface area in m²
    bbox: BoundingBox
        The w/e/s/n bounding box to fill with the grid

    Returns
    =======
    grid: numpy.array, shape=(M, N, 2)
        one rectangular sub-grid of the hexagonal grid
    mid_grid:
        the complementary rectangular sub-grid of the hexagonal grid

    """

    # Just in case, deal with bounding boxes crossing the date line. This may
    # not be fine for the plotter, I don't know, but the logic should be sound,
    # I hope.
    # TODO: Write a regression test for that case.

    if bbox.e < bbox.w:
        east = bbox.e + 360
    else:
        east = bbox.e
    bbox_centre = Point(
        (east + bbox.w)/2,
        (bbox.n + bbox.s)/2)
    hexagon_side, grid_point_distance = hexagon(area)

    points = [bbox_centre]

    # Neighbors east and west: Direct tiling
    while points[-1].longitude < east:
        next = GEODESIC.direct(points[-1], 90, grid_point_distance)
        points.append(Point(*numpy.array(next)[0, :2]))
    while points[0].longitude > bbox.w:
        next = GEODESIC.direct(points[0], 270, grid_point_distance)
        points.insert(0, Point(*numpy.array(next)[0, :2]))

    grid = numpy.array([points])
    while (grid[0, :, 1] < bbox.n).any():
        next = GEODESIC.direct(grid[0], 0, 3 * hexagon_side)
        grid = numpy.vstack((
            [numpy.array(next)[:, :2]],
            grid))
    while (grid[-1, :, 1] > bbox.s).any():
        next = GEODESIC.direct(grid[-1], 180, 3 * hexagon_side)
        grid = numpy.vstack((
            grid,
            [numpy.array(next)[:, :2]]))

    mid_grid = (grid[:-1, :-1, :] + grid[1:, 1:, :]) / 2

    return numpy.array((grid[:-1, :-1, :], mid_grid))


# Define continents

area = 450_000_000 #m²

def is_land(xy):
   return LAND.contains(sgeom.Point(*xy))

def coordinates_to_index(points, resolution=2 * 60):
    """Convert long,lat coordinate pairs into indices in a TIF

    Convert a [..., 2] ndarray, or a pair of coordinates, into the matching
    grid indices of a Mercator projection pixel map with a given resolution in
    pixels per degree.

    Paramaters
    ==========
    points: ndarray-like, shape=(..., 2)
        An array of longitude. latitude pairs to be converted to grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray(int), shape=(..., 2)
        An integer array of grid indices

    """
    points = numpy.asarray(points)
    return numpy.stack(
        (numpy.round((-points[..., 1] + 90) * resolution).astype(int),
         numpy.round((points[..., 0] + 180) * resolution).astype(int)),
        -1)

def index_to_coordinates(indices, resolution=2 * 60):
    """Convert grid indices into long,lat coordinate pairs

    Convert a (ndarray of) grid indices of a Mercator projection pixel map with
    a given resolution in pixels per degree into geocoordinate pairs (longitude,
    latitude).

    Paramaters
    ==========
    indices: ndarray(int), shape=(..., 2)
        An integer array of grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray, shape=(..., 2)
        An array of longitude. latitude pairs

    """
    indices = numpy.asarray(indices)
    return numpy.stack(
        (indices[..., 1] / resolution - 180,
         90 - indices[..., 0] / resolution),
        -1)

# TODO: Write tests for middle, random, each corner, forwards and backwards.


class Grid():
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


    def __init__(self, m, i, j):
        self.m = m
        self.ij = i, j
        self.population = 0
        self.popcap = self.model.population_capacity(
            self.point) * area / 1000000 # area is in m², popcap is per km²
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


class PopulationCapModel:
    alpha = 10 ** -8.07
    beta = 2.64

    @property
    def precipitation_tif(self):
        try:
            return self._tif
        except AttributeError:
            self._tif = tifffile.imread("../worldclim/wc2.0_bio_30s_12.tif").clip(0)
            return self._tif

    def population_capacity(self, point):
        """Calculate the pop cap of a cell given its precipitation

        Return the carrying capacity K of a hexagonal cell of area AREA with
        the given mean yearly precipitation measured in mm

        In Gavin et al. 2017, the function used is

        K = α * P ** β

        with eg. α = 10 ** -8.07, β = 2.64 [Note: Their supplementary material
        lists α=-8.07, but that is a nonsensical number. The actual plot line
        does not exactly correspond to α=10^-8.07, but more closely to
        α=10^-7.96, but that suggests that this is at least close to the
        described behaviour.]

        Parameters
        ==========
        precipitation: float
            The cell's mean yearly precipitation P, in mm

        Returns
        =======
        float
            The cell's carrying capacity, in individuals/km^2
        """
        return self.alpha * self.precipitation(point) ** self.beta

    def precipitation(self, point):
        index = tuple(coordinates_to_index(point))
        return self.precipitation_tif[index]




namerica = BoundingBox(
    w=-178.8,
    s=6.6,
    e=-49.0,
    n=83.3)
australia = BoundingBox(
    112.8708,
    153.7392,
    -43.8615,
    -9.6712)
americas = BoundingBox(
    e=-34.535395,
    s=-56.028198,
    w=-168.571541,
    n=74.526716
)


# Start the simulation
if __name__ == "__main__":
    from language import DifferenceSpreadAndRoundLanguage as Language

    run = hex(numpy.random.randint(4096))

    continent = americas

    # Generate a hexagonal grid over the continent
    class LGrid(Grid):
        grid = hexagonal_earth_grid(continent, area)
        all_gridcells = {}
        model = PopulationCapModel()

    land = (
        numpy.apply_along_axis(is_land, axis=2, arr=LGrid.grid[0]),
        numpy.apply_along_axis(is_land, axis=2, arr=LGrid.grid[1]))

    # Find a random starting cell that is capable of supporting at least one individual
    g = LGrid.random_cell()
    while LGrid.gridcell(*g).popcap < 1:
        g = LGrid.random_cell()

    # Create a population in that starting grid cell
    l = Language(1, LGrid.gridcell(*g))
    generation = 0
    while True:
        generation += 1
        while True:
            try:
                l.grow()
            except StopIteration as s:
                print([c.population for c in l.cells])
                print(l.popcap)
                print(s)
                break
        filled = {g for g in LGrid.all_gridcells.values()
                if g.language is not None}
        expansion_zone = {i for g in filled
                        for i in g.neighbors()
                        if i.popcap > 0} - filled
        if not expansion_zone:
            break
        free = list(expansion_zone)
        l = Language(l.id + 1, free[numpy.random.randint(len(free))])


        # Plot the results
        plt.gcf().set_size_inches(15, 15)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines("50m")
        ax.set_extent(continent)

        colors = {}
        ax.add_collection(plot_hex_grid(
            lambda cell: colors.setdefault(
                cell.language, numpy.random.random(size=3))
            if cell.language else (0, 0, 0, 0),
            LGrid.all_gridcells))
        plt.savefig("output_{:}_{:08}.png".format(run, generation),
                    dpi=300)
        plt.close()

    def l_id(l_or_none):
        if l_or_none is None:
            return None
        if l_or_none.language is None:
            return None
        return l_or_none.language.id

    json.dump(
        [
            [[l_id(all_gridcells.get((0, i, j)))
              for i in range(LGrid.grid[0].shape[0])]
             for j in range(LGrid.grid[0].shape[1])],
            [[l_id(all_gridcells.get((1, i, j)))
              for i in range(LGrid.grid[1].shape[0])]
             for j in range(LGrid.grid[1].shape[1])]
        ],
        open("output_{:}.json".format(run), "w"))

    plt.gcf().set_size_inches(15, 15)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines("50m")
    ax.set_extent(continent)

    colors = {}
    ax.add_collection(plot_hex_grid(
        lambda cell: colors.setdefault(
            cell.language, numpy.random.random(size=3))
        if cell.language else (0, 0, 0, 0),
        LGrid.all_gridcells))
    plt.savefig("output_{:}_final.png".format(run),
                dpi=300)
    plt.close()

