import sys
import argparse
from pathlib import Path

import numpy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot simulation results")
parser.add_argument("logfile", type=argparse.FileType("r"))
parser.add_argument("output_dir", type=Path, default="plots/", nargs="?")
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--max-cdif", type=int, default=20,
                    help="Maxiumum cultural difference to be plotted")
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)

def bitvec_to_color(i: int):
    r = (bin(i & 0b000000000000111111).count("1") - 1) / 5
    g = (bin(i & 0b000000111111000000).count("1") - 1) / 5
    b = (bin(i & 0b111111000000000000).count("1") - 1) / 5
    r, g, b = min(max(r, 0), 1), min(max(g, 0), 1), min(max(b, 0), 1)
    return 1-r, 1-g, b


stem = Path(args.logfile.name).stem
all_q = numpy.zeros((51, args.max_cdif))
m = 0
n = 0
pop = []
try:
    for l, line in enumerate(args.logfile):
        if line.startswith("POPULATION: ["):
            m += 1
            if m < args.start:
                continue
            content = eval(line[len("POPULATION: "):])
            population = [(x, y, sum(p.values()))
                          for x, y, p in content]
            colors = [bitvec_to_color(max(p, key=p.get))
                      for x, y, p in content]
            q = numpy.array(population).T
            q[:2] *= 180 / numpy.pi
            pop.append(q[2].sum())
            q[2] /= 4  # Scale population to pixels
            plt.scatter(*q[:2], q[2], c=colors, alpha=0.8, linewidths=0.0)
            plt.xlim(-170, -30)
            plt.ylim(-60, 90)
            plt.gcf().set_size_inches((12, 16))
            plt.savefig(args.output_dir / f"disp{m:08d}-{stem:}.png")
            plt.close()
        if line.startswith("GD_CD: {"):
            n += 1
            if n < args.start:
                continue
            print("gd_cd {:d}".format(n))
            line = eval(line[len("GD_CD: "):])
            q = numpy.zeros_like(all_q)
            for i in range(q.shape[0]):
                for j in range(q.shape[1]):
                    q[i, j] = line.get(i, {}).get(j, 0)
            q /= q.max(1)[:, None] + 1
            plt.imshow(q.T)
            all_q += q
            plt.gcf().set_size_inches((12, 9))
            plt.savefig(args.output_dir / f"corr{n:08d}-{stem:}.png")
            plt.close()
finally:
    plt.imshow(all_q.T)
    plt.gcf().set_size_inches((12, 9))
    plt.savefig(args.output_dir / f"corr-{stem:}.png")
    plt.close()

    plt.plot(pop)
    plt.savefig(args.output_dir / f"pop-{stem:}.png")
    plt.close()
