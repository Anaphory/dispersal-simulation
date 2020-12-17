import sys
import json
import argparse
import itertools
import typing as t
from pathlib import Path

import numpy
from numpy import pi as π, cos
from cartopy.geodesic import Geodesic
import datetime
import matplotlib.pyplot as plt

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

# graph = {
#     (min(i, j), max(i, j)): d
#     for (i, j, d) in json.load(args.statefile)["p"]["dispersal_graph"]["edges"]
# }


def bitvec_to_color(i: int):
    r = (bin(i & 0b000000000000111111).count("1") - 1) / 5
    g = (bin(i & 0b000000111111000000).count("1") - 1) / 5
    b = (bin(i & 0b111111000000000000).count("1") - 1) / 5
    r, g, b = min(max(r, 0), 1), min(max(g, 0), 1), min(max(b, 0), 1)
    return 1 - r, 1 - g, b


for logfile in args.logfile:
    stem = Path(logfile.name).stem
    if stem == "<stdin>":
        stem = hex(int(datetime.datetime.now().timestamp()))
    all_q = numpy.zeros((51, args.max_cdif))
    m = 0
    n = 0
    pop = []
    try:
        for l, line in enumerate(logfile):
            if line.startswith("POPULATION: ["):
                m += 1
                if m < args.start:
                    continue
                content = eval(line[len("POPULATION: ") :])
                population = [(x, y, sum(p.values())) for x, y, p in content]
                colors = [bitvec_to_color(max(p, key=p.get)) for x, y, p in content]
                q = numpy.array(population).T
                pop.append(q[2].sum())
                q[2] /= 4  # Scale population to pixels
                plt.scatter(*q[:2], q[2], c=colors, alpha=0.8, linewidths=0.0)
                plt.xlim(*args.xlim)
                plt.ylim(*args.ylim)
                plt.gcf().set_size_inches((12, 16))
                plt.savefig(args.output_dir / f"disp{m:08d}-{stem:}.png")
                print(args.output_dir / f"disp{m:08d}-{stem:}.png")
                plt.close()

        scatter: t.Dict[t.Tuple[int, int], int] = {}
        for ((x1, y1, p1), (x2, y2, p2)) in itertools.combinations_with_replacement(content, 2):
            # if abs(y1-y2) > 2.7:
            #     continue
            # if abs(x1-x2) > 2.7/cos(π*(y1+y2)/2):
            #     continue
            flat_dist = int(g.inverse((x1, y1), (x2, y2))[0, 0] / 1000 + 0.5)
            # if flat_dist > 300000:
            #     continue
            for c1, n1 in p1.items():
                for c2, n2 in p2.items():
                    cult_dist = bin(c1 ^ c2).count("1")
                    scatter.setdefault((flat_dist, cult_dist), 0)
                    scatter[flat_dist, cult_dist] += n1 * n2
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
    finally:
        plt.plot(pop)
        plt.savefig(args.output_dir / f"pop-{stem:}.png")
        plt.close()
