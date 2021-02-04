import sys
import json
import argparse
import itertools
import typing as t
from pathlib import Path
from tqdm import tqdm

import numpy
from numpy import pi as π, cos
from cartopy.geodesic import Geodesic
import datetime
import matplotlib.pyplot as plt

import osm

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
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)


def bitvec_to_color(i: int):
    r = (bin(i & 0b000000000000111111).count("1") - 1) / 5
    g = (bin(i & 0b000000111111000000).count("1") - 1) / 5
    b = (bin(i & 0b111111000000000000).count("1") - 1) / 5
    r, g, b = min(max(r, 0), 1), min(max(g, 0), 1), min(max(b, 0), 1)
    return 1 - r, 1 - g, b

season_length_in_years = 1/6.

regions=list(osm.areas)
def ccontains(polygon):
    cache = {}
    def is_in(location):
        try:
            return cache[location]
        except KeyError:
            cache[location] = osm.contains(
                polygon,
                location)
            return cache[location]
    return is_in
cf =  [ccontains(osm.shapes[a]) for a in regions]


def compute_contained_population(population, containment_functions=cf):
    pop = [0 for _ in containment_functions]

    for x, y, p, c in population:
        for i, inside in enumerate(containment_functions):
            if inside((x, y)):
                pop[i] += p
    return pop


distances = {}


for logfile in args.logfile:
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
        if line.startswith("Parameters"):
            i0 = line.index("season_length_in_years")
            before = line.index(" ", i0)
            after = line.index(" ", before + 1)
            season_length_in_years = float(line[before:after -1])
            continue
        elif line.startswith("t: "):
            timestamp = int(line[3:]) * season_length_in_years
        elif line.startswith("POPULATION: ["):
            m += 1
            if m < args.start:
                continue
            content = [(x, y, n, c)
                    for x, y, p in eval(line[len("POPULATION: ") :])
                    for c, n in p.items()]
            # smallest to be plotted last:
            content.sort(
                key=lambda xync: xync[2],
                reverse=True)
            xs, ys, ns, cs = zip(*content)
            ns = numpy.array(ns) # Scale population to pixels
            plt.scatter(xs, ys, ns / 4, c=[bitvec_to_color(c) for c in cs], alpha=0.5, linewidths=0.0)
            plt.xlim(*args.xlim)
            plt.ylim(*args.ylim)
            plt.gcf().set_size_inches((12, 16))
            plt.savefig(args.output_dir / f"disp{m:08d}-{stem:}.png")
            print(args.output_dir / f"disp{m:08d}-{stem:}.png")
            plt.close()

            ts.append(timestamp)
            subpops.append(compute_contained_population(content))
            pop.append(ns.sum())

    plt.gcf().set_size_inches((16, 9))
    plt.gca().set_yscale("log")
    plt.plot(ts, pop, label="Total population")
    for subpop, r in zip(zip(*subpops), regions):
        plt.plot(ts, subpop, label=r)
    plt.legend()
    plt.xlim(0, 2000)
    plt.savefig(args.output_dir / f"pop-{stem:}.png")
    plt.close()

    plt.gcf().set_size_inches((12, 9))
    plt.gca().set_yscale("log")
    plt.plot(ts, pop, label="Total population")
    for subpop, r in zip(zip(*subpops), regions):
        plt.plot(ts, subpop, label=r)
    plt.legend()
    plt.savefig(args.output_dir / f"pop-{stem:}.png")
    plt.close()

    scatter: t.Dict[t.Tuple[int, int], int] = {}
    for ((x1, y1, p1, c1), (x2, y2, p2, c2)) in tqdm(itertools.combinations_with_replacement(content, 2), total = len(content) * (len(content)  + 1) // 2):
        try:
            flat_dist = distances[(x1, x2, y1, y2)]
        except KeyError:
            flat_dist = int(g.inverse((x1, y1), (x2, y2))[0, 0] / 1000 + 0.5)
            distances[(x1, x2, y1, y2)] = flat_dist
        cult_dist = bin(c1 ^ c2).count("1")
        scatter.setdefault((flat_dist, cult_dist), 0)
        scatter[flat_dist, cult_dist] += p1 * p2
    try:
        x, y, count = zip(*[(f, c, n) for (f, c), n in scatter.items()])
        count = numpy.array(count, dtype=float)
        count *= 400 / count.max()
    except ValueError:
        continue
    plt.scatter(x, y, s=count)
    plt.gcf().set_size_inches((12, 9))
    plt.savefig(args.output_dir / f"corr{m:08d}-{stem:}.png")
    plt.close()
