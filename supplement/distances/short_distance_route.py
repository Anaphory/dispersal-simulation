import itertools
from itertools import count
import typing as t
from heapq import heappush as push, heappop as pop

from tqdm import tqdm
from sqlalchemy import func
import numpy
from numpy import pi, cos, inf
from sqlalchemy import select

import rasterio
import cartopy.io.shapereader as shpreader
import shapely.wkb

import cartopy.crs as ccrs
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap


from raster_data import ecoregion_tile, boundingbox_from_tile, Tile, RowCol
from database import db
from earth import GEODESIC
from main import load_distances

DATABASE, TABLES = db("sqlite:///migration-network.sqlite")


def distances_from_focus(
    source: RowCol,
    destination: RowCol,
    trafo: rasterio.Affine,
    distance_by_direction: t.Dict[RowCol, numpy.array],
    pred: t.Optional = None,
    max_time: float = 3.0 * 3600,
) -> numpy.array:
    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    d = distance_by_direction[0, 1]
    dist: numpy.array = numpy.full((d.shape[0], d.shape[1] + 1), numpy.inf, dtype=float)
    seen: t.Dict[t.Tuple[int, int], float] = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples.
    # Speeds faster than 2m/s are really rare, even on rivers, so use geodesic
    # distance with a speed of 2m/s as heuristic
    heuristic = (
        GEODESIC.inverse(trafo * source[::-1], trafo * destination[::-1])[0, 0] / 4
    )
    fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []
    push(fringe, (heuristic, 0, next(c), source))

    def moore_neighbors(
        r0: int, c0: int
    ) -> t.Iterable[t.Tuple[t.Tuple[int, int], float]]:
        for (r, c), d in distance_by_direction.items():
            r1, c1 = r0 + r, c0 + c
            if 0 <= r1 < d.shape[0] and 0 <= c1 < d.shape[1]:
                yield (r1, c1), min(max_time, d[r1, c1])

    while fringe:
        (_, d, _, spot) = pop(fringe)
        if dist[spot] < d:
            continue  # already searched this node.
        dist[spot] = d
        if spot == destination:
            break

        for u, cost in moore_neighbors(*spot):
            vu_dist = dist[spot] + cost
            if numpy.isfinite(dist[u]) and vu_dist < dist[u]:
                print(
                        f"Contradictory paths found: Already reached u={u} at distance {dist[u]} via {pred[u]} [{dist[pred[u]]}], but now finding a shorter connection {vu_dist} via {spot} [{dist[spot]}]. Do you have negative weights, or a bad heuristic?"
                    )
            if u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                heuristic = (
                    GEODESIC.inverse(trafo * (u[::-1]), trafo * (destination[::-1]))[0, 0]
                    / 4
                )
                push(fringe, (vu_dist + heuristic, vu_dist, next(c), u))
                if pred is not None:
                    pred[u] = spot
    return dist


def pixel(lon, lat) -> RowCol:
    col, row = ~trafo * (lon, lat)
    return int(row + 0.5), int(col + 0.5)


distances, trafo = load_distances(("N", 30, "W", 90))

nashville: RowCol = pixel(-86.9759, 36.0346)
natchez: RowCol = pixel(-91.36892, 31.54543)

# nashville: RowCol = (3872, 1248)
# natchez: RowCol = (3872, 983)
print(nashville, natchez)

# Plotting
ax = plt.axes()

print("Nashville to Natchez…")
pred = {}
d = distances_from_focus(nashville, natchez, trafo, distances, pred=pred)

backtrack = natchez
rows = [natchez[0]]
cols = [natchez[1]]
print("Backtrack…")
while pred.get(backtrack):
    backtrack = pred[backtrack]
    rows.append(backtrack[0])
    cols.append(backtrack[1])

ax.plot(cols, rows, zorder=30, linewidth=4, c="blue")

print("Natchez to Nashville…")
rpred = {}
rd = distances_from_focus(natchez, nashville, trafo, distances, pred=rpred)
print(rd[nashville])
rbacktrack = nashville
rrows = [nashville[0]]
rcols = [nashville[1]]
print("Backtrack…")
while rpred.get(rbacktrack):
    rbacktrack = rpred[rbacktrack]
    rrows.append(rbacktrack[0])
    rcols.append(rbacktrack[1])

minrow = min(*rows, *rrows)
maxrow = max(*rows, *rrows)
mincol = min(*cols, *rcols)
maxcol = max(*cols, *rcols)

margin = 10

print(minrow, maxrow, mincol, maxcol)
# minrow, maxrow, mincol, maxcol = 3834, 3910, 983, 1260

ax.set_xlim(mincol - 0.5 - margin, maxcol + 0.5 + margin)
ax.set_ylim(maxrow + 0.5 + margin, minrow - 0.5 - margin)

print("Set up plot…")
ax.plot(rcols, rrows, zorder=31, linewidth=3, c="red")

isolatedness = numpy.min(
    numpy.stack(
        [
            d[
                minrow - margin : maxrow + 1 + margin,
                mincol - margin : maxcol + 1 + margin,
            ]
            for d in distances.values()
        ]
    ),
    axis=0,
)
plt.imshow(
    isolatedness,
    extent=(
        mincol - 0.5 - margin,
        maxcol + 0.5 + margin,
        maxrow + 0.5 + margin,
        minrow - 0.5 - margin,
    ),
    origin="upper",
    zorder=10,
)

print("Show…")
plt.show()
