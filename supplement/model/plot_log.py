import json
import numpy
import igraph
import itertools
import matplotlib
import cartopy.crs as ccrs
from matplotlib import pyplot as plt

import gavin2017processbased
from util import hexagon_coords, neighbors
from dispersal import cultural_distance


collection = None

def connect_agents(family_locations, G=None):
    if G is None:
        G = igraph.Graph()

    xyz = {}
    edges = []
    for location1, families1 in family_locations:
        if not families1:
            continue
        coord = coordinates[tuple(location1)]
        for (name, size, _) in families1:
            G.add_vertices(name)
            lon, lat = numpy.random.normal(0, 0.04, 2) + coord
            xyz[len(xyz)] = (lon, lat,  size / 5)
        for location2, families2 in family_locations:
            if not families2:
                continue
            if abs(location1[1] - location2[1]) > 6:
                continue
            max_lat = max(abs(location1[1]), abs(location2[1]))
            if ((location1[0] - location2[0]) * numpy.cos(max_lat / 180 * numpy.pi)) ** 2 + (location1[1] - location2[1]) ** 2 > 6 ** 2:
                continue
            edges.append(connect_similar_agents(families1, families2))
            if location2 == location1:
                break
    G.add_edges(itertools.chain(*edges))
    return G, xyz

def connect_similar_agents(families1, families2):
    for (id1, size1, culture1) in families1:
        for (id2, size2, culture2) in families2:
            if id1 == id2:
                break
            if cultural_distance(culture1, culture2) > 6:
                continue
            yield id1, id2

def plot(family_locations, resources, maximum=3650000000.0, hexes=[], color_schema={}) -> None:
    plt.gcf().set_size_inches(30, 30)
    ax = plt.axes(projection=ccrs.PlateCarree())

    if resources is not None:
        if not hexes:
            for index in numpy.ndindex(*coordinates.shape[:-1]):
                hexes.append(numpy.array(hexagon_coords(index, neighbors(index), coordinates)))
            global collection
        collection = matplotlib.collections.PolyCollection(hexes)

        cmap = plt.get_cmap("viridis")
        vmax = max(resources)
        values = [cmap(v / vmax)[:-1] + (0.2,) if v > 0 else [0, 0, 0, 0] for v in resources]
        collection.set_facecolor(values)
        ax.add_collection(collection)

    ax.coastlines("50m")
    ax.set_extent(gavin2017processbased.americas)

    G, xyz = connect_agents(family_locations)

    new_color_schema = {}
    # for community in G.community_label_propagation():
    # for community in G.community_infomap():
    for community in G.community_fastgreedy().as_clustering():
        pos = numpy.array([(xyz[n][0], xyz[n][1]) for n in community])
        s = [xyz[n][2] for n in community]
        mean = numpy.mean(pos, axis=0)
        closest = numpy.inf
        for location, popsize in color_schema:
            dist = numpy.linalg.norm(mean - location) + abs(numpy.log((popsize + 1) / (len(community) + 1))) / numpy.log(1.04)
            if dist > 20:
                continue
            if dist < closest:
                closest = dist
                reference = location, popsize
        if numpy.isfinite(closest):
            color = color_schema.pop(reference)
        else:
            color = numpy.random.random(3)

        color

        new_color_schema[tuple(mean), len(community)] = color

        plt.scatter(*pos.T,
                s=s,
                c=[color],
                edgecolors=None)

    for key in list(color_schema):
        color_schema.pop(key)
    color_schema.update(new_color_schema)


def plot_(s):
    plot([(family.effective_size, location, family.culture)
          for location, families in s.families.items()
          for family in families],
         [s.grid.patches[mij].resources
          for mij in numpy.ndindex(*coordinates.shape[:-1])])

    plt.show()

def plot_last_line(filename):
    with open(filename, "r") as f:
        for line in f:
            pass
    properties = json.loads(line)
    # DAFUQ? Why even are there 00 bytes??
    families = properties["Families"].values()
    hexes = properties["Resources"]
    plot(families, hexes)
    plt.show()

def plot_series(f, template="dispersal-{:07d}.png", limit=None, show_resources=True):
    for l, line in enumerate(f):
        if limit is not None and l not in limit:
            continue
        print(l)
        properties = json.loads(line)
        families = properties["Families"]
        hexes = properties.get("Resources", None)
        plot(families, show_resources if hexes else None)
        plt.savefig(template.format(l))
        plt.close()


def plot_alaska_population(filename):
    pop = []
    cache = {}
    def is_in(location):
        try:
            return cache[location]
        except KeyError:
            cache[location] = osm.contains(osm.shapes["Alaska"], coordinates[location])
            return cache[location]
    with open(filename, "r") as f:
        for l, line in enumerate(f):
            pop.append(0)
            try:
                data = json.loads(line)
            except json.JSONDecodeError:
                break
            for effective_size, location, culture in data["Families"].values():
                if is_in(tuple(location)):
                    pop[l] += effective_size
    plt.plot(pop)
    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=argparse.FileType("r"), default="log")
    parser.add_argument("--template", default="plots/dispersal-{:07d}.pdf")
    parser.add_argument("--area", type=int, default=450_000_000)
    args = parser.parse_args()
    coordinates: numpy.ndarray = gavin2017processbased.hexagonal_earth_grid(
        gavin2017processbased.americas,
        args.area)
    plot_series(args.log, template=args.template, show_resources=False)
