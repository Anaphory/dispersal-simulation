import sys
import json
import argparse
import itertools
import typing as t
from pathlib import Path
from tqdm import tqdm
import re

import numpy
from numpy import pi as π, cos
from cartopy.geodesic import Geodesic
import datetime
import matplotlib.pyplot as plt
from matplotlib import cm
from colorspacious import cspace_converter

import osm
import networkx

g = Geodesic()

parser = argparse.ArgumentParser(description="Plot simulation results")
# parser.add_argument(
#     "statefile",
#     type=argparse.FileType("r"),
#     help="JSON file containing the state, and thus the graph",
# )
parser.add_argument("logfile", type=argparse.FileType("r"), nargs="+")
parser.add_argument("output_dir", type=Path, default="plots/", nargs="?")
parser.add_argument("--start", type=int, default=0)
parser.add_argument(
    "--max-cdif",
    type=int,
    default=20,
    help="Maxiumum cultural difference to be plotted",
)
parser.add_argument(
    "--xlim", type=lambda x: [int(i) for i in x.split(":")], default="-170:-30"
)
parser.add_argument(
    "--ylim", type=lambda x: [int(i) for i in x.split(":")], default="-60:90"
)
parser.add_argument(
    "--start-year",
    "-t",
    default=1400,
    type=int,
    help="First year where the whole continent is settled, so it works as a starting point",
)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)


def bitvec_to_color(i: int):
    r = (bin(i & 0b000000000000111111).count("1") - 1) / 5
    g = (bin(i & 0b000000111111000000).count("1") - 1) / 5
    b = (bin(i & 0b111111000000000000).count("1") - 1) / 5
    r, g, b = min(max(r, 0), 1), min(max(g, 0), 1), min(max(b, 0), 1)
    return 1 - r, 1 - g, b


season_length_in_years = 1 / 6.0

regions = list(osm.areas)


def ccontains(polygon):
    cache = {}

    def is_in(location):
        try:
            return cache[location]
        except KeyError:
            cache[location] = osm.contains(polygon, location)
            return cache[location]

    return is_in


cf = [ccontains(osm.shapes[a]) for a in regions]


def compute_contained_population(population, containment_functions=cf):
    pop = [0 for _ in containment_functions]

    for (x, y), p in population:
        for i, inside in enumerate(containment_functions):
            if inside((x, y)):
                pop[i] += p
    return pop


CUTOFF = 30

for logfile in args.logfile:
    popcaps = {}
    actual_pops = {}
    stem = Path(logfile.name).stem
    if stem == "<stdin>":
        stem = hex(int(datetime.datetime.now().timestamp()))
    all_q = numpy.zeros((51, args.max_cdif))
    m = 0
    n = 0
    ts = []
    pop = []
    subpops = []
    for l, line in enumerate(logfile):
        if line.startswith(" Parameters"):
            (param,) = re.findall("season_length_in_years: *([0-9.]+)", line)
            season_length_in_years = float(param)
            edges = eval(re.findall(r"edges: ([(), 0-9]*),", line)[0])
            nodes = eval(
                re.findall(
                    r"node weights: *(\{([0-9]+: *\(.*\{[^}]*\} *\)[, ]*)+\})", line
                )[0][0]
            )
            nodes_from_coords = {(x, y): n for n, (h, x, y, p) in nodes.items()}
            G = networkx.Graph()
            G.add_edges_from(edges)
            print(G)
            continue
        elif line.startswith("Resources"):
            try:
                x, y, p = re.match(
                    "Resources .*? \((-?[0-9.]+), (-?[0-9.]+)\), .*popcap ([0-9.]+) .*",
                    line,
                ).groups()
            except AttributeError:
                continue
            loc = float(x), float(y)
            popcaps[loc] = float(p)
            actual_pops[loc] = []
        elif line.startswith("t: "):
            timestamp = int(line[3:]) * season_length_in_years
        elif line.startswith("POPULATION: ["):
            m += 1
            if m < args.start:
                continue
            content = [
                (x, y, n, c)
                for x, y, p in eval(line[len("POPULATION: ") :])
                for c, n in p.items()
            ]
            # smallest to be plotted last:
            content.sort(key=lambda xync: xync[2], reverse=True)
            xs, ys, ns, cs = zip(*content)
            ns = numpy.array(ns)  # Scale population to pixels
            plt.scatter(
                xs,
                ys,
                ns / 4,
                c=[bitvec_to_color(c) for c in cs],
                alpha=0.5,
                linewidths=0.0,
            )
            plt.xlim(*args.xlim)
            plt.ylim(*args.ylim)
            plt.gcf().set_size_inches((12, 16))
            plt.savefig(args.output_dir / f"disp{m:08d}-{stem:}.png")
            print(args.output_dir / f"disp{m:08d}-{stem:}.png")
            plt.close()

            ts.append(timestamp)
            subpops.append(
                compute_contained_population(((x, y), p) for x, y, p, c in content)
            )
            pop.append(ns.sum())
            if timestamp >= args.start_year:
                sns = {
                    (x, y): sum(n for c, n in p.items())
                    for x, y, p in eval(line[len("POPULATION: ") :])
                }
                for loc, nn in sns.items():
                    try:
                        actual_pops[loc].append(nn)
                    except KeyError:
                        print(f"Did not find {loc}.")

    caps = compute_contained_population(popcaps.items())

    print("Population development in the first 2000 years…")
    l = 0
    plt.gcf().set_size_inches((16, 9))
    plt.gca().set_yscale("log")
    for subpop, c, r in zip(zip(*subpops), caps, regions):
        plt.plot(ts, subpop, label=r, c=cm.tab20(l))
        plt.plot([0, 2000], [c, c], c=cm.tab20(l + 1), linestyle="dotted")
        l += 2
    plt.plot(ts, pop, label="Total population", c=cm.tab20(l))
    plt.legend()
    plt.xlim(0, 2000)
    plt.savefig(args.output_dir / f"pop2k-{stem:}.png")
    plt.close()

    print("Population development over the whole run…")
    l = 0
    plt.gcf().set_size_inches((12, 9))
    plt.gca().set_yscale("log")
    for subpop, c, r in zip(zip(*subpops), caps, regions):
        plt.plot(ts, subpop, label=r, c=cm.tab20(l))
        plt.plot([0, max(ts)], [c, c], c=cm.tab20(l + 1), linestyle="dotted")
        l += 2
    plt.plot(ts, pop, label="Total population", c=cm.tab20(l))
    plt.legend()
    plt.savefig(args.output_dir / f"pop-{stem:}.png")
    plt.close()

    print("Population caps and actual populations in each spot…")
    mean_actual_pops = {loc: (sum(pop)/len(pop) if pop else 0) for loc, pop in actual_pops.items()}
    plt.gcf().set_size_inches((16, 9))
    plt.scatter(popcaps.values(), mean_actual_pops.values(), s=4, alpha=0.03, c='k')
    # plt.plot(range(int(max(popcaps.values())+0.5)), mean_aggregated_popsize_bands, c="0.5")
    plt.xlim(0, 600)
    plt.ylim(0, 600)
    plt.plot((0, 600), (0, 600))
    plt.savefig(args.output_dir / f"popcap-{stem:}.png")
    plt.close()

    print("Population caps and actual populations in each spot and its immediate neigbors…")
    mean_actual_pops_n0 = {
        nodes_from_coords[loc]: mean_pop for loc, mean_pop in mean_actual_pops.items()
    }
    popcaps_n0 = {
        nodes_from_coords[loc]: mean_pop for loc, mean_pop in popcaps.items()
    }
    mean_actual_pops_n1 = mean_actual_pops_n0.copy()
    popcaps_n1 = popcaps_n0.copy()
    mean_actual_pops_n2 = mean_actual_pops_n0.copy()
    popcaps_n2 = popcaps_n0.copy()

    for node in G:
        try:
            n2 = {node}
            for neighbor in G[node]:
                mean_actual_pops_n1[node] = mean_actual_pops_n1[node] + mean_actual_pops_n0.get(neighbor, 0)
                popcaps_n1[node] = popcaps_n1[node] + popcaps_n0.get(neighbor, 0)
                n2.add(neighbor)
                mean_actual_pops_n2[node] = mean_actual_pops_n2[node] + mean_actual_pops_n0.get(neighbor, 0)
                popcaps_n2[node] = popcaps_n2[node] + popcaps_n0.get(neighbor, 0)
                for neighbor2 in G[neighbor]:
                    if neighbor2 not in n2:
                        n2.add(neighbor2)
                        mean_actual_pops_n2[node] = mean_actual_pops_n2[node] + mean_actual_pops_n0.get(neighbor2, 0)
                        popcaps_n2[node] = popcaps_n2[node] + popcaps_n0.get(neighbor2, 0)
            mean_actual_pops_n2[node] /= len(n2)
            popcaps_n2[node] /= len(n2)
        except KeyError:
            continue

    print("Computing color representation of points…")
    lab = cspace_converter("CIELab", "sRGB1")
    colors = numpy.clip([
        lab((40, n[1]-10, n[0]+79)) for n in popcaps.keys()
    ], 0, 1)

    print("Plotting Color Scheme…")
    plt.gcf().set_size_inches((12, 16))
    plt.scatter(*zip(*popcaps.keys()), s=10, c=colors,
        alpha=0.6,
        linewidths=0.0,
    )
    plt.xlim(*args.xlim)
    plt.ylim(*args.ylim)
    plt.savefig(args.output_dir / f"colorschema-{stem}.png")
    plt.close()

    print("Plotting intended vs. opserved population…")
    plt.gcf().set_size_inches((16, 9))
    plt.scatter(popcaps.values(), mean_actual_pops.values(), s=2, alpha=0.4, linewidths=0.0, c=colors)
    plt.xlim(0, 400)
    plt.ylim(0, 400)
    plt.plot((0, 400), (0, 400))
    plt.savefig(args.output_dir / f"popcap-{stem:}.png")
    plt.close()

    plt.gcf().set_size_inches((16, 9))
    plt.scatter(popcaps_n1.values(), mean_actual_pops_n1.values(), s=2, alpha=0.4, linewidths=0.0, c=colors)
    plt.xlim(0, 2000)
    plt.ylim(0, 2000)
    plt.plot((0, 2000), (0, 2000))
    plt.savefig(args.output_dir / f"popcap_n1-{stem:}.png")
    plt.close()

    plt.gcf().set_size_inches((16, 9))
    plt.scatter(popcaps_n2.values(), mean_actual_pops_n2.values(), s=1, alpha=0.8, linewidths=0.0, c=colors)
    plt.xlim(0, 250)
    plt.ylim(0, 250)
    plt.plot((0, 500), (0, 500))
    plt.savefig(args.output_dir / f"popcap_n2-{stem:}.png")
    plt.close()

    print("Computing geographic distance vs. cultural distance…")
    start = None
    scatter: t.Dict[t.Tuple[int, int], int] = {}
    scattere: t.Dict[t.Tuple[int, int], int] = {}
    for i, ((x1, y1, p1, c1), (x2, y2, p2, c2)) in tqdm(enumerate(
        itertools.combinations(content, 2)), total=len(content) * (len(content) - 1) // 2
    ):
        if i%23 != 0:
            continue
        if start != nodes_from_coords[x1, y1]:
            start = nodes_from_coords[x1, y1]
            edge_dists = networkx.single_source_shortest_path_length(
                G, start, cutoff=CUTOFF
            )
        flat_dist = int(g.inverse((x1, y1), (x2, y2))[0, 0] / 1000 + 0.5)
        edge_dist = edge_dists.get(nodes_from_coords[x2, y2], CUTOFF + 2)
        cult_dist = bin(c1 ^ c2).count("1")
        scatter.setdefault((flat_dist, cult_dist), 0)
        scatter[flat_dist, cult_dist] += p1 * p2
        scattere.setdefault((edge_dist, cult_dist), 0)
        scattere[edge_dist, cult_dist] += p1 * p2

    print("Geodesic distance vs. cultural distance…")
    try:
        x, y, count = zip(*[(f, c, n) for (f, c), n in scatter.items()])
        count = numpy.array(count, dtype=float)
        count *= 25 / count.max()
    except ValueError:
        continue
    plt.scatter(x, y, s=count)
    plt.gcf().set_size_inches((12, 9))
    plt.savefig(args.output_dir / f"corr{m:08d}-{stem:}.png")
    plt.close()

    print("Edge count vs. cultural distance…")
    try:
        x, y, count = zip(*[(f, c, n) for (f, c), n in scattere.items()])
        x = numpy.array(x)
        count = numpy.array(count, dtype=float)
        for x0 in set(x):
            count[x == x0] *= 100 / count[x == x0].max()
    except ValueError:
        continue
    plt.scatter(x, y, s=count)
    plt.gcf().set_size_inches((12, 9))
    plt.savefig(args.output_dir / f"corre{m:08d}-{stem:}.png")
    plt.close()
