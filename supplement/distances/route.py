import itertools
from itertools import count
import typing as t
from heapq import heappush as push, heappop as pop

from sqlalchemy import func
import numpy
from numpy import pi, cos, inf
from sqlalchemy import select

import rasterio
import cartopy.geodesic as geodesic
import cartopy.io.shapereader as shpreader
import shapely.wkb

from database import db
from earth import LAND, BBOX, GEODESIC

DATABASE, TABLES = db("sqlite:///migration-network-backup.sqlite")


def distances_from_focus(
    source,
    destination,
    filter_sources: t.Optional[t.Set[str]] = None,
    pred: t.Optional = None,
) -> numpy.array:
    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    seen: t.Dict[t.Tuple[int, int], float] = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples
    fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []
    push(fringe, (0, next(c), source))

    if filter_sources:
        filter = [TABLES["edges"].c.source.in_(filter_sources)]
    else:
        filter = []

    dist = t.DefaultDict(lambda: inf)

    while fringe:
        (d, _, spot) = pop(fringe)
        if dist.get(spot) is not None:
            continue  # already found this node, faster, through another path
        dist[spot] = d
        if spot == destination:
            break

        for u, cost, source in itertools.chain(
            DATABASE.execute(
                select(
                    TABLES["edges"].c.node2,
                    TABLES["edges"].c.travel_time,
                    TABLES["edges"].c.source,
                ).where(TABLES["edges"].c.node1 == spot, *filter)
            ),
            # TODO: Remove reverse direction!
            DATABASE.execute(
                select(
                    TABLES["edges"].c.node1,
                    TABLES["edges"].c.travel_time,
                    TABLES["edges"].c.source,
                ).where(TABLES["edges"].c.node2 == spot, *filter)
            ),
        ):
            vu_dist = dist[spot] + cost
            if u in dist and vu_dist < dist[u]:
                raise ValueError(
                    "Contradictory paths found. Do you have negative weights?"
                )
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if pred is not None:
                    pred[u] = spot, source
    return dist


def find_node(lon, lat, rivers=False):
    """Find the closest node near given geocoordinates.

    If rivers=True, include river nodes.
    """

    # Find the seven nearest nodes, in terms of raw coordinates, using an approximation
    query = (
        select(
            TABLES["nodes"].c.node_id,
            TABLES["nodes"].c.longitude,
            TABLES["nodes"].c.latitude,
        )
        .where(
            TABLES["nodes"].c.longitude != None,
            TABLES["nodes"].c.latitude != None,
        )
        .order_by(
            func.abs(TABLES["nodes"].c.latitude - lat)
            + func.abs(TABLES["nodes"].c.longitude - lon) * cos(lat * pi / 180)
        )
        .limit(7)
    )

    if not rivers:
        query = query.where(TABLES["nodes"].c.node_id > 100000000)

    min_d = numpy.inf
    for node, nlon, nlat in DATABASE.execute(query):
        d = GEODESIC.inverse((lon, lat), (nlon, nlat))[0, 0]
        if d < min_d:
            best = node
            min_d = d

    return best


nashville = find_node(-86.9759, 36.0346)
natchez = find_node(-91.36892, 31.54543)

# Plotting
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.feature import ShapelyFeature
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap


proj = ccrs.PlateCarree()
ax = plt.axes(projection=proj)

natchez_trace_parkway = shapely.wkb.load(open("./natchez_trace.wkb", "rb"))
ax.add_geometries(
    natchez_trace_parkway,
    ccrs.PlateCarree(),
    facecolor="green",
    alpha=0.5,
    edgecolor=[0, 0.5, 0, 1],
    zorder=32,
)

pred = {nashville: (None, "self-loop")}
d = distances_from_focus(nashville, natchez, pred=pred, filter_sources={"grid"})
backtrack = natchez
lons = [-91.36892]
lats = [31.54543]
sources = []
while pred.get(backtrack):
    backtrack, source = pred[backtrack]
    backtrack = pred[backtrack]
    lon, lat = DATABASE.execute(
        select(
            TABLES["nodes"].c.longitude,
            TABLES["nodes"].c.latitude,
        ).where(TABLES["nodes"].c.node_id == backtrack)
    ).fetchone()
    lons.append(lon)
    lats.append(lat)
    sources.append(source)
bbox = (min(lons), max(lons), min(lats), max(lats))
ax.plot(lons, lats, zorder=30)

rpred = {natchez: (None, "self-loop")}
rd = distances_from_focus(natchez, nashville, pred=rpred)
rbacktrack = nashville
rlons = [-86.9759]
rlats = [36.0346]
rsources = []
while rpred.get(rbacktrack):
    backtrack, source = rpred[rbacktrack]
    lon, lat = DATABASE.execute(
        select(
            TABLES["nodes"].c.longitude,
            TABLES["nodes"].c.latitude,
        ).where(TABLES["nodes"].c.node_id == rbacktrack)
    ).fetchone()
    rlons.append(lon)
    rlats.append(lat)
    rsources.append(source)

bbox = (
    min(bbox[0], *rlons),
    max(bbox[1], *rlons),
    min(bbox[2], *rlats),
    max(bbox[3], *rlats),
)
ax.plot(rlons, rlats, zorder=30)

ax.set_extent(bbox)  # x0, x1, y0, y1

random = ListedColormap(
    [(1, 1, 1, 0)] + [numpy.random.random(3) for _ in range(2 ** 16 - 1)], name="random"
)

for tile in itertools.product(
    ["N", "S"], [10, 30, 50, 70], ["E", "W"], [30, 60, 90, 120, 150, 180]
):
    fname_v = "voronoi-{:s}{:d}{:s}{:d}.tif".format(*tile)
    try:
        data = rasterio.open(fname_v, "r").read(1)[500:-500, 500:-500]
    except rasterio.errors.RasterioIOError:
        continue
    west = (-1 if tile[2] == "W" else 1) * tile[3]
    south = (-1 if tile[0] == "S" else 1) * tile[1]
    v = ax.imshow(
        data,
        transform=proj,
        extent=(west, west + 30, south, south + 20),
        interpolation="nearest",
        cmap=random,
        vmin=0,
        vmax=2 ** 16 - 1,
    )
    v.set_zorder(1)
    # TODO: we tend to run out of memory otherwise, but at least test!

# shape_feature = ShapelyFeature(LAND, ccrs.PlateCarree(), facecolor="none", edgecolor="blue", lw=1)
ax.add_geometries(
    shpreader.Reader(
        shpreader.natural_earth(resolution="10m", category="physical", name="lakes")
    ).geometries(),
    ccrs.PlateCarree(),
    facecolor="white",
    alpha="0.8",
    edgecolor="none",
    zorder=2,
)

ax.add_geometries(
    shpreader.Reader(
        shpreader.natural_earth(resolution="10m", category="physical", name="ocean")
    ).geometries(),
    ccrs.PlateCarree(),
    facecolor="white",
    edgecolor="none",
    zorder=2,
)

hexagon_id, x, y, hx, hy, coast = zip(
    *DATABASE.execute(
        select(
            TABLES["nodes"].c.node_id,
            TABLES["nodes"].c.longitude,
            TABLES["nodes"].c.latitude,
            TABLES["nodes"].c.h3longitude,
            TABLES["nodes"].c.h3latitude,
            TABLES["nodes"].c.coastal,
        ).where(TABLES["nodes"].c.node_id > 100000000)
    ).fetchall()
)


ax.scatter(
    hx,
    hy,
    marker=".",
    edgecolors="none",
    facecolors=[{False: [0, 0, 0, 1], True: [0.5, 0.5, 0.5, 0.9]}[c] for c in coast],
    s=9,
    zorder=3,
)
ax.scatter(
    x,
    y,
    marker=".",
    facecolors="none",
    edgecolors=[{False: [0, 0, 0, 1], True: [0.5, 0.5, 0.5, 0.9]}[c] for c in coast],
    zorder=3,
)
node_id, x, y, coast = zip(
    *DATABASE.execute(
        select(
            TABLES["nodes"].c.node_id,
            TABLES["nodes"].c.longitude,
            TABLES["nodes"].c.latitude,
            TABLES["nodes"].c.coastal,
        ).where(TABLES["nodes"].c.node_id <= 100000000)
    ).fetchall(),
)
ax.scatter(
    x,
    y,
    c="teal",
    marker="+",
    s=25,
    zorder=4,
)

plt.show()
