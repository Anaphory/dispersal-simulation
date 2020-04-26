import json
import itertools

import numpy
import igraph
import matplotlib
import cartopy.crs as ccrs
from matplotlib import pyplot as plt

import dispersal_model.osm as osm
import dispersal_model.hexgrid as hexgrid
from dispersal_model.dispersal import cultural_distance


collection = None


def connect_agents(family_locations, G=None):
    if G is None:
        G = igraph.Graph()

    xyz = {}
    edges = []
    for location1, families1 in family_locations:
        if not families1:
            continue
        coord = hexgrid.geo_coordinates(location1)
        for (name, size, _) in families1:
            G.add_vertices(name)
            lon, lat = numpy.random.normal(0, 0.01, 2) + coord
            xyz[len(xyz)] = (lon, lat,  size)
        for location2, families2 in family_locations:
            if not families2:
                continue
            if hexgrid.hex_distance(location1, location2) > 8:
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


def plot(family_locations,
         resources, maximum=3650000000.0,
         hexes=[],
         color_schema={}) -> None:
    plt.gcf().set_size_inches(30, 30)
    ax = plt.axes(projection=ccrs.Sinusoidal(-98))

    if resources is not None:
        raise NotImplementedError

    ax.coastlines("50m")
    ax.set_extent((hexgrid.AMERICAS.w,
                   hexgrid.AMERICAS.e,
                   hexgrid.AMERICAS.s,
                   hexgrid.AMERICAS.n))

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
            dist = numpy.linalg.norm(mean - location) + abs(
                numpy.log((popsize + 1) / (len(community) + 1))) / numpy.log(
                    1.04)
            if dist > 20:
                continue
            if dist < closest:
                closest = dist
                reference = location, popsize
        if numpy.isfinite(closest):
            color = color_schema.pop(reference)
        else:
            color = numpy.random.random(3)

        new_color_schema[tuple(mean), len(community)] = color

        ax.scatter(*pos.T,
                   s=s,
                   c=[color],
                   alpha=0.2,
                   edgecolors=None,
                   transform=ccrs.Geodetic())

    for key in list(color_schema):
        color_schema.pop(key)
    color_schema.update(new_color_schema)


def plot_series(f, template="dispersal-{:07d}.png",
                limit=None, show_resources=True):
    for l, line in enumerate(f):
        if limit is not None and l not in limit:
            continue
        print(l)
        properties = json.loads(line)
        families = properties["Families"]
        hexes = properties.get("Resources", None)
        plot(families, hexes if show_resources else None)
        plt.savefig(template.format(l))
        plt.close()


from h3 import h3


def plot_alaska_population(filename):
    pop = []
    cache = {}

    def is_in(location):
        try:
            return cache[location]
        except KeyError:
            cache[location] = osm.contains(
                osm.shapes["Alaska"],
                hexgrid.geo_coordinates(location))
            return cache[location]

    with open(filename, "r") as f:
        for l, line in enumerate(f):
            pop.append(0)
            try:
                data = json.loads(line)
            except json.JSONDecodeError:
                break
            for location, families in data["Families"]:
                if is_in(location):
                    for descendence, effective_size, culture in families:
                        pop[l] += effective_size
    plt.plot(pop)
    plt.show()


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=argparse.FileType("r"), default="log")
    parser.add_argument("--template", default="plots/dispersal-{:07d}.pdf")
    parser.add_argument("--area", type=int, default=450_000_000)
    parser.add_argument("--limit", type=eval, default=None)
    parser.add_argument("--resources", action="store_true", default=False)
    args = parser.parse_args()
    plot_series(args.log,
                template=args.template,
                limit=args.limit,
                show_resources=args.resources)


if __name__ == "__main__":
    main()
