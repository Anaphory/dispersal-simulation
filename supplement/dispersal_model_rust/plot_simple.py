import re
import sys
import argparse
from pathlib import Path

import numpy
import matplotlib.pyplot as plt
from matplotlib.cm import Set3

parser = argparse.ArgumentParser(description="Plot simulation results")
parser.add_argument("logfile", type=argparse.FileType("r"))
parser.add_argument("output_dir", type=Path, default="plots/", nargs="?")
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--max-cdif", type=int, default=20,
                    help="Maxiumum cultural difference to be plotted")
parser.add_argument("--xlim", type=lambda x: [int(i) for i in x.split(":")], default="-170:-30")
parser.add_argument("--ylim", type=lambda x: [int(i) for i in x.split(":")], default="-60:90")
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)


stem = Path(args.logfile.name).stem
all_q = numpy.zeros((51, args.max_cdif))
m = 0
n = 0
pop = []
target_pop = []
for l, line in enumerate(args.logfile):
    if line.startswith("Resources at"):
        popcap = int(re.search("popcap (\\d+)", line).group(1))
        target_pop.append(popcap)
        pop.append([])
    elif line.startswith("POPULATION: ["):
        m += 1
        if m < args.start:
            continue
        content = eval(line[len("POPULATION: "):])
        population = {-y: sum(p.values()) for x, y, p in content}
        for i in range(len(pop)):
            pop[i].append(population.get(i, 0))

plt.gca().set_yscale("log")
for i, ns in enumerate(pop):
    plt.scatter(0, numpy.mean(ns[len(ns)//2:])/target_pop[i], color=Set3(i/len(pop)), s=120)
    plt.scatter(numpy.arange(len(ns)),  numpy.array(ns)/target_pop[i], color=Set3(i/len(pop)), label=str(target_pop[i]))
plt.legend()
plt.savefig(args.output_dir / f"pop-{stem:}.png")
plt.show()
plt.close()
