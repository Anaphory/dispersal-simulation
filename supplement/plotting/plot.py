#!/home/cluster/gkaipi/.pyenv/shims/python
# -*- encoding: utf-8 -*-
from ast import literal_eval
import sys
import json
import shutil
import argparse
import itertools
import typing as t
from pathlib import Path
from tqdm import tqdm
import re
import scipy.signal
import collections

import numpy
from numpy import pi, cos
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
    "--last-only",
    action="store_true",
    default=False,
    help="investigate only the last two years",
)
parser.add_argument(
    "--start-year",
    "-t",
    default=1400,
    type=int,
    help="First year where the whole continent is settled, so it works as a starting point",
)
parser.add_argument(
    "--graph",
    action="store_true",
    default=False,
    help="Parse the graph structure from the log file and use it to compute neighbor statistics",
)
parser.add_argument(
    "--all-steps",
    action="store_true",
    default=False,
    help="Generate all intermediate dispersal plots",
)
parser.add_argument(
    "--gd-cd",
    action="store_true",
    default=False,
    help="Generate georaphical distance/cultural distance plots",
)
args = parser.parse_args()
try:
    args.output_dir.mkdir(parents=True)
except OSError:
    pass


def bitvec_to_color(i, top=[0, 0, 0]):
    rgb = (
        bin(i)[-1:0:-3].count("1"),
        bin(i)[-2:0:-3].count("1"),
        bin(i)[-3:0:-3].count("1"),
    )
    top[0] = max(rgb[0] + 1, top[0])
    top[1] = max(rgb[1] + 1, top[1])
    top[2] = max(rgb[2], top[2])
    return [rgb[i] / top[i] for i in range(3)]


def cultural_distance(c1, c2):
    return bin(c1 ^ c2).count("1")


season_length_in_years = 1 / 6.0

regions = list(osm.areas)
tdf = regions.index("Tierra del Fuego (Isla Grande)")


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


def count_cultures(culture_set, _cache={}):
    culture_set = frozenset(culture_set)
    try:
        return _cache[culture_set]
    except KeyError:
        pass
    culture_groups = []
    for c1 in culture_set:
        groups = [
            (g, group)
            for g, group in enumerate(culture_groups)
            if any(
                cultural_distance(c1, c2) < params["cooperation_threshold"]
                for c2 in group
            )
        ]
        if len(groups) == 0:
            culture_groups.append({c1})
        else:
            groups[0][1].add(c1)
            for g, group in groups[:0:-1]:
                groups[0][1].update(group)
                del culture_groups[g]
    _cache[culture_set] = len(culture_groups)
    return _cache[culture_set]


def compute_contained_population(population, last_reference, containment_functions=cf):
    pop = [0 for _ in containment_functions]
    culture = [set() for _ in containment_functions]

    for (x, y), (p, c) in population:
        for i, inside in enumerate(containment_functions):
            if inside((x, y)):
                pop[i] += p
                culture[i].add(c)

    n_cultures_now = []
    times = []
    for i, ((time, cs_then), cs_now) in enumerate(zip(last_reference, culture)):
        n_now = count_cultures(cs_now)
        n_cultures_now.append(n_now)
        n_then = count_cultures(cs_then)
        n_both = count_cultures(cs_then | cs_now)
        if n_both == n_now + n_then:
            # There is no continuous culture from then to now
            last_reference[i] = [0, frozenset(cs_now)]
        else:
            last_reference[i][0] += 1
        times.append(last_reference[i][0])

    return pop, n_cultures_now, times


correction = {}
def find_nearby(everywhere, loc):
    s = numpy.inf
    for x, y in everywhere:
        if abs(x-loc[0]) + abs(y-loc[1]) < s:
            nearby = x, y
            s = abs(x-loc[0]) + abs(y-loc[1]) 
    return nearby

def plot_content(content, ts, pop, subpops, cultures_by_location, persistence_by_location, last_references, plot=True):
    # smallest to be plotted last:
    content.sort(key=lambda xync: xync[2], reverse=True)
    xs, ys, ns, cs = zip(*content)
    ns = numpy.array(ns)  # Scale population to pixels
    if plot:
        plt.text(
            0.08,
            0.99,
            "{:7d}".format(int(timestamp)),
            fontsize="large",
            ha="right",
            va="top",
            transform=plt.gca().transAxes,
        )
        plt.scatter(
            xs,
            ys,
            numpy.maximum(ns / 400, 4),
            c=[bitvec_to_color(c) for c in cs],
            alpha=0.5,
            linewidths=0.0,
        )
        plt.xlim(*args.xlim)
        plt.ylim(*args.ylim)
        plt.gcf().set_size_inches((12, 16))
        plt.savefig(
            str(args.output_dir / "disp{m:08d}-{stem:}.png".format(m=m, stem=stem)),
            bbox_inches="tight",
        )
        plt.close()
    print(
        str(
            (
                args.output_dir / "disp{m:08d}-{stem:}.png".format(m=m, stem=stem)
            ).absolute()
        )
    )

    ts.append(timestamp)
    pop.append(ns.sum())
    populations, cultures, times = compute_contained_population((
        ((x, y), (p, c)) for x, y, p, c in content), last_references
    )
    subpops.append(populations)
    cultures_by_location.append(cultures)
    persistence_by_location.append(times)

    sns = {
        (x, y): sum(n for x, y, n, c in p)
        for (x, y), p in itertools.groupby(content, key=lambda x: (x[0], x[1]))
    }
    for loc, nn in sns.items():
        try:
            actual_pops[correction.get(loc, loc)].append(nn)
        except KeyError:
            correction[loc] = find_nearby(actual_pops, loc)
            print("Did not find {loc}.".format(loc=loc))


CUTOFF = 30
bitvec = re.compile('"\[([01]*)\]"')

for logfile in args.logfile:
    print("Starting {:} at {:}…".format(logfile, datetime.datetime.now()), file=sys.stderr)
    params = {
        "resource_recovery_per_season": 0.10,
        "culture_mutation_rate": 6e-3,
        "culture_dimensionality": 20,
        "cooperation_threshold": 6,
        "maximum_resources_one_adult_can_harvest": 0.25,
        "evidence_needed": 0.1,
        "payoff_std": 0.1,
        "minimum_adaptation": 0.5,
        "fight_deadliness": (2 ** 31),
        "enemy_discount": 0.5 ** (0.2),
        "season_length_in_years": 1.0 / 6.0,
        "warfare": "true",
        "end": "???",
        "last": 0,
        "mean_pop": 0.0,
    }

    bitvec_to_color.__defaults__[0][0] = 1
    bitvec_to_color.__defaults__[0][1] = 1
    bitvec_to_color.__defaults__[0][2] = 1
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
    cultures_by_location = []
    persistence_by_location = []
    last_reference = [[0, frozenset()] for _ in regions]
    family_path = collections.defaultdict(list)
    for l, line in enumerate(logfile):
        if line.startswith(" Parameters"):
            igraph = line.index("dispersal_graph:")
            bracket_level = 0
            for i, s in enumerate(line[igraph:], igraph):
                if s == "{":
                    bracket_level += 1
                elif s == "}":
                    bracket_level -= 1
                    if bracket_level == 0:
                        break
            end_graph = i + 1
            param_string = (
                line[len(" Parameters { ") : igraph] + line[end_graph : -len(" }\n")]
            )
            for param_value in param_string.split(","):
                param_value = param_value.strip()
                if not param_value:
                    continue
                param, value = param_value.split(":")
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                params[param] = value
            season_length_in_years = params["season_length_in_years"]
            if args.last_only:
                recent_contents = [
                    None for _ in range(int(2 / season_length_in_years) + 1)
                ]

            edges = literal_eval(
                re.findall(r"edges: ([(), 0-9]*),", line[igraph:end_graph])[0]
            )
            nodes = literal_eval(
                re.findall(
                    r"node weights: *(\{([0-9]+: *\(.*\{[^}]*\} *\)[, ]*)+\})",
                    line[igraph:end_graph],
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
            if args.last_only:
                recent_contents.pop()
                recent_contents.insert(0, line[len("POPULATION: ") :])
                continue
            try:
                line = bitvec.sub("0b\\1", line)
                content = [
                    (x, y, n, c)
                    for x, y, p in literal_eval(line[len("POPULATION: ") :])
                    for c, n in p.items()
                ]
            except SyntaxError:
                print("Could not parse line", line)
                continue
            if timestamp < args.start:
                continue
            plot_content(
                content, ts, pop, subpops, cultures_by_location, persistence_by_location, last_reference, plot=args.all_steps
            )
        elif line.startswith("F:"):
            match = re.search(r'escendence: "([^"]*)".*ocation: NodeIndex\(([0-9]*)\),.*culture: ([01]*)', line)
            if match:
                name, node, culture = match.group(1), match.group(2), match.group(3)
                _, x, y, _ = nodes[int(node)]
                family_path[name].append((x, y))
            else:
                print(line)
    if line.startswith("Ended"):
        params["end"] = "Ended"
    elif line.startswith("Died out"):
        params["end"] = "Died out"
    elif line.startswith("slurmstepd:"):
        params["end"] = "Timeout"
    else:
        params["end"] = "???"

    if args.last_only:
        m -= len(recent_contents)
        for line in recent_contents:
            if not line:
                continue
            m += 1
            content = [
                (x, y, n, c) for x, y, p in literal_eval(line) for c, n in p.items()
            ]
            plot_content(content, ts, pop, subpops, cultures_by_location, persistence_by_location, last_reference, True)
    elif not args.all_steps:
        plot_content(content, ts, pop, subpops, cultures_by_location,  persistence_by_location, last_reference, True)
    shutil.copyfile(
        (args.output_dir / "disp{m:08d}-{stem:}.png".format(m=m, stem=stem)),
        (args.output_dir / "disp-last-{stem:}.png".format(stem=stem)),
    )
    print(args.output_dir / "disp-last-{stem:}.png".format(stem=stem))

    caps, _, _ = compute_contained_population(
        [((x, y), (p, 1)) for (x, y), p in popcaps.items()], [[0, frozenset()] for _ in last_reference]
    )

    params["mean_pop"] = sum(pop) / len(pop)
    params["last"] = timestamp
    for subpop, culturecount, persistence, c, r in zip(
            zip(*subpops), zip(*cultures_by_location), zip(*persistence_by_location), caps, regions
    ):
        params[f"{r}_relative"] = numpy.mean(subpop) / c
        params[f"{r}_persistence"] = numpy.mean(numpy.diff(scipy.signal.find_peaks(persistence, width=1)[0]))
        params[f"{r}_cultures"] = numpy.median(
            [c for c, n in zip(culturecount, subpop) if n > 0]
        )
        try:
            params[f"{r}_arrival"] = [x > 0 for x in subpop].index(True)
        except ValueError:
            params[f"{r}_arrival"] = numpy.nan

    summaries = Path("runs_overview.tsv")
    if not summaries.exists():
        with summaries.open("w") as paramfile:
            print("Path", *params.keys(), sep="\t", file=paramfile)
    with summaries.open("a") as paramfile:
        print(Path(logfile.name).absolute(), *params.values(), sep="\t", file=paramfile)
        print(Path(logfile.name).absolute(), *params.values(), sep="\t")

    print("Migration of some families…")
    for family in family_path.values():
        jittered = (numpy.random.random(size=(len(family), 2))-0.5) * [0.01, 0.02] + family
        plt.plot(*zip(*jittered),  linewidth=1)
        plt.scatter(*zip(*jittered),  alpha=0.2)
    plt.xlim(*args.xlim)
    plt.ylim(*args.ylim)
    plt.gcf().set_size_inches((12, 16))
    plt.savefig(
        str(args.output_dir / "migration-{stem:}.png".format(m=m, stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print("Population development in the first 2000 years…")
    l = 0
    plt.gcf().set_size_inches((22, 12))
    plt.gca().set_yscale("log")
    for subpop, c, r in zip(zip(*subpops), caps, regions):
        plt.plot(ts, subpop, label=r, c=cm.Set3(l))
        plt.plot([0, 2000], [c, c], c=cm.Set3(l), alpha=0.5, linestyle="dotted")
        l += 1
    plt.plot(ts, pop, label="Total population", c="k")
    plt.legend()
    plt.xlim(0, 2000)
    plt.savefig(
        str(args.output_dir / "pop2k-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print("Population development in the last 2000 years…")
    l = 0
    plt.gcf().set_size_inches((22, 12))
    plt.gca().set_yscale("log")
    for subpop, c, r in zip(zip(*subpops), caps, regions):
        plt.plot(ts, subpop, label=r, c=cm.Set3(l))
        plt.plot(
            [max(ts) - 2000, max(ts)],
            [c, c],
            c=cm.Set3(l),
            alpha=0.5,
            linestyle="dotted",
        )
        l += 1
    plt.plot(ts, pop, label="Total population", c="k")
    plt.legend()
    plt.xlim(max(ts) - 2000, max(ts))
    plt.savefig(
        str(args.output_dir / "popl2k-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print("Population development over the whole run…")
    l = 0
    plt.gcf().set_size_inches((22, 12))
    plt.gca().set_yscale("log")
    for subpop, c, r in zip(zip(*subpops), caps, regions):
        plt.plot(ts, subpop, label=r, c=cm.Set3(l))
        plt.plot([0, max(ts)], [c, c], c=cm.Set3(l), linestyle="dotted")
        l += 1
    plt.plot(ts, pop, label="Total population", c="k")
    plt.legend()
    plt.savefig(
        str(args.output_dir / "pop-{stem:}.png".format(stem=stem)), bbox_inches="tight"
    )
    plt.close()

    print("Culture counts over the whole run…")
    l = 0
    plt.gcf().set_size_inches((22, 12))
    filter_width = min(19, len(cultures_by_location) // 5 * 2 + 1)
    for subpop, c, r in zip(zip(*cultures_by_location), caps, regions):
        plt.plot(ts, scipy.signal.medfilt(subpop, filter_width), label=r, c=cm.Set3(l))
        l += 1
    plt.legend()
    plt.savefig(
        str(args.output_dir / "cultcount-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print("Culture counts over the last 2000 years…")
    l = 0
    plt.gcf().set_size_inches((22, 12))
    filter_width = min(11, len(cultures_by_location) // 5 * 2 + 1)
    for subpop, c, r in zip(zip(*cultures_by_location), caps, regions):
        plt.plot(ts, scipy.signal.medfilt(subpop, filter_width), label=r, c=cm.Set3(l))
        l += 1
    plt.legend()
    plt.xlim(max(ts) - 2000, max(ts))
    plt.savefig(
        str(args.output_dir / "cultcountl2k-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print("Culture persistence over the whole run…")
    l = 0
    plt.gcf().set_size_inches((22, 12))
    filter_width = min(19, len(cultures_by_location) // 5 * 2 + 1)
    for subpop, c, r in zip(zip(*persistence_by_location), caps, regions):
        plt.plot(ts, subpop, label=r, c=cm.Set3(l))
        l += 1
    plt.legend()
    plt.savefig(
        str(args.output_dir / "cultper-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print("Culture persistence over the last 2000 years…")
    l = 0
    plt.gcf().set_size_inches((22, 12))
    filter_width = min(11, len(cultures_by_location) // 5 * 2 + 1)
    for subpop, c, r in zip(zip(*persistence_by_location), caps, regions):
        plt.plot(ts, subpop, label=r, c=cm.Set3(l))
        l += 1
    plt.legend()
    plt.xlim(max(ts) - 2000, max(ts))
    plt.savefig(
        str(args.output_dir / "cultper2k-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()


    print("Population caps and actual populations in each spot…")
    mean_actual_pops = {
        loc: (sum(pop) / len(pop) if pop else 0) for loc, pop in actual_pops.items()
    }
    plt.gcf().set_size_inches((16, 9))
    plt.scatter(popcaps.values(), mean_actual_pops.values(), s=4, alpha=0.03, c="k")
    # plt.plot(range(int(max(popcaps.values())+0.5)), mean_aggregated_popsize_bands, c="0.5")
    plt.xlim(0, 600)
    plt.ylim(0, 600)
    plt.plot((0, 600), (0, 600))
    plt.savefig(
        str(args.output_dir / "popcap-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print(
        "Population caps and actual populations in each spot and its immediate neigbors…"
    )
    mean_actual_pops_n0 = {
        nodes_from_coords[loc]: mean_pop for loc, mean_pop in mean_actual_pops.items()
    }
    popcaps_n0 = {nodes_from_coords[loc]: mean_pop for loc, mean_pop in popcaps.items()}

    print("Computing color representation of points…")
    lab = cspace_converter("CIELab", "sRGB1")
    colors = numpy.clip([lab((40, n[1] - 10, n[0] + 79)) for n in popcaps.keys()], 0, 1)

    print("Plotting Color Scheme…")
    plt.gcf().set_size_inches((12, 16))
    x, y = zip(*popcaps.keys())
    plt.scatter(x, y, s=10, c=colors, alpha=0.6, linewidths=0.0)
    plt.xlim(*args.xlim)
    plt.ylim(*args.ylim)
    plt.savefig(
        str(args.output_dir / "colorschema-{stem}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()

    print("Plotting intended vs. observed population…")
    plt.gcf().set_size_inches((16, 9))
    plt.scatter(
        popcaps.values(),
        mean_actual_pops.values(),
        s=2,
        alpha=0.4,
        linewidths=0.0,
        c=colors,
    )
    plt.xlim(0, 400)
    plt.ylim(0, 400)
    plt.plot((0, 400), (0, 400))
    plt.savefig(
        str(args.output_dir / "popcap-{stem:}.png".format(stem=stem)),
        bbox_inches="tight",
    )
    plt.close()


    if args.graph:
        mean_actual_pops_n1 = mean_actual_pops_n0.copy()
        popcaps_n1 = popcaps_n0.copy()
        mean_actual_pops_n2 = mean_actual_pops_n0.copy()
        popcaps_n2 = popcaps_n0.copy()

        G = networkx.Graph()
        edges = literal_eval(
            re.findall(r"edges: ([(), 0-9]*),", graph_string)[0]
        )
        nodes = literal_eval(
            re.findall(
                r"node weights: *(\{([0-9]+: *\(.*\{[^}]*\} *\)[, ]*)+\})",
                graph_string,
            )[0][0]
        )
        nodes_from_coords = {(x, y): n for n, (h, x, y, p) in nodes.items()}
        G.add_edges_from(edges)
        print(G)

        for node in G:
            try:
                n2 = {node}
                for neighbor in G[node]:
                    mean_actual_pops_n1[node] = mean_actual_pops_n1[
                        node
                    ] + mean_actual_pops_n0.get(neighbor, 0)
                    popcaps_n1[node] = popcaps_n1[node] + popcaps_n0.get(neighbor, 0)
                    n2.add(neighbor)
                    mean_actual_pops_n2[node] = mean_actual_pops_n2[
                        node
                    ] + mean_actual_pops_n0.get(neighbor, 0)
                    popcaps_n2[node] = popcaps_n2[node] + popcaps_n0.get(neighbor, 0)
                    for neighbor2 in G[neighbor]:
                        if neighbor2 not in n2:
                            n2.add(neighbor2)
                            mean_actual_pops_n2[node] = mean_actual_pops_n2[
                                node
                            ] + mean_actual_pops_n0.get(neighbor2, 0)
                            popcaps_n2[node] = popcaps_n2[node] + popcaps_n0.get(
                                neighbor2, 0
                            )
                mean_actual_pops_n2[node] /= len(n2)
                popcaps_n2[node] /= len(n2)
            except KeyError:
                continue

        plt.gcf().set_size_inches((16, 9))
        plt.scatter(
            popcaps_n1.values(),
            mean_actual_pops_n1.values(),
            s=2,
            alpha=0.4,
            linewidths=0.0,
            c=colors,
        )
        plt.xlim(0, 2000)
        plt.ylim(0, 2000)
        plt.plot((0, 2000), (0, 2000))
        plt.savefig(
            str(args.output_dir / "popcap_n1-{stem:}.png".format(stem=stem)),
            bbox_inches="tight",
        )
        plt.close()

        plt.gcf().set_size_inches((16, 9))
        plt.scatter(
            popcaps_n2.values(),
            mean_actual_pops_n2.values(),
            s=1,
            alpha=0.8,
            linewidths=0.0,
            c=colors,
        )
        plt.xlim(0, 250)
        plt.ylim(0, 250)
        plt.plot((0, 500), (0, 500))
        plt.savefig(
            str(args.output_dir / "popcap_n2-{stem:}.png".format(stem=stem)),
            bbox_inches="tight",
        )
        plt.close()

        if not args.gd_cd:
            continue

        print("Computing geographic distance vs. cultural distance…")
        start = None
        scatter = {}
        scattere = {}
        for i, ((x1, y1, p1, c1), (x2, y2, p2, c2)) in tqdm(
            enumerate(itertools.combinations(content, 2)),
            total=len(content) * (len(content) - 1) // 2,
        ):
            if i % 23 != 0:
                continue
            if start != nodes_from_coords[x1, y1]:
                start = nodes_from_coords[x1, y1]
                edge_dists = networkx.single_source_shortest_path_length(
                    G, start, cutoff=CUTOFF
                )
            flat_dist = int(g.inverse((x1, y1), (x2, y2))[0, 0] / 1000 + 0.5)
            edge_dist = edge_dists.get(nodes_from_coords[x2, y2], CUTOFF + 2)
            cult_dist = cultural_distance(c1, c2)
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
        plt.savefig(
            str(args.output_dir / "corr{m:08d}-{stem:}.png".format(m=m, stem=stem)),
            bbox_inches="tight",
        )
        plt.close()

        print("Distribution of neighbor families")
        y, count = zip(*[(c, n) for (f, c), n in scattere.items() if f == 1])
        plt.bar(y, count)
        plt.gcf().set_size_inches((12, 9))
        plt.savefig(
            str(args.output_dir / "neighb{m:08d}-{stem:}.png".format(m=m, stem=stem)),
            bbox_inches="tight",
        )
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
        plt.savefig(
            str(args.output_dir / "corre{m:08d}-{stem:}.png".format(m=m, stem=stem)),
            bbox_inches="tight",
        )
        plt.close()
        print("Finished {:} at {:}.".format(logfile, datetime.datetime.now()), file=sys.stderr)
