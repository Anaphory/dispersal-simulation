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
from matplotlib import cm

from database import db
from earth import GEODESIC

DATABASE, TABLES = db("sqlite:///migration-network.sqlite")


def distances_from_focus(
    source,
    destination,
    filter_sources: t.Optional[t.Set[str]] = None,
    filter_nodes: bool = False,
    pred: t.Optional = None,
) -> numpy.array:
    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html,
    # extended to be an A* implementation
    seen: t.Dict[t.Tuple[int, int], float] = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples
    fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []
    source_lonlat = DATABASE.execute(
        select(
            [
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
            ]
        ).where(TABLES["nodes"].c.node_id == source)
    ).fetchone()
    dest_lonlat = DATABASE.execute(
        select(
            [
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
            ]
        ).where(TABLES["nodes"].c.node_id == destination)
    ).fetchone()
    # Speeds faster than 3m/s are really rare, even on rivers, so use geodesic
    # distance with a speed of 3m/s as heuristic
    # heuristic = GEODESIC.inverse(source_lonlat, dest_lonlat)[0, 0] / 3
    heuristic = 0
    push(fringe, (heuristic, 0, next(c), source))

    if filter_sources:
        filter = TABLES["edges"].c.source.in_(filter_sources)
        if filter_nodes:
            filter &= TABLES["nodes"].c.node_id > 100000000
    else:
        if filter_nodes:
            filter = TABLES["nodes"].c.node_id > 100000000
        else:
            filter = True

    dist = t.DefaultDict(lambda: inf)

    while fringe:
        (_, d, _, spot) = pop(fringe)
        if dist.get(spot, numpy.inf) < d:
            continue  # already found this node, faster, through another path
        dist[spot] = d
        if spot == destination:
            break

        for u, cost, source, lon, lat in DATABASE.execute(
            select(
                [
                    TABLES["edges"].c.node2,
                    TABLES["edges"].c.travel_time,
                    TABLES["edges"].c.source,
                    TABLES["nodes"].c.longitude,
                    TABLES["nodes"].c.latitude,
                ]
            )
            .select_from(
                TABLES["edges"].join(
                    TABLES["nodes"],
                    onclause=TABLES["edges"].c.node2 == TABLES["nodes"].c.node_id,
                )
            )
            .where(
                (
                    (
                        (TABLES["edges"].c.node1 == spot)
                        & (TABLES["edges"].c.node2 != spot)
                    )
                    if filter is True
                    else filter
                )
                & (TABLES["edges"].c.node1 == spot)
                & (TABLES["edges"].c.node2 != spot)
            )
        ):
            vu_dist = dist[spot] + cost
            if u in dist and vu_dist < dist[u]:
                print(
                    f"Contradictory paths found: Already reached u={u} at distance {dist[u]} via {pred[u]} [{dist[pred[u]]}], but now finding a shorter connection {vu_dist} via {spot} [{dist[spot]}]. Do you have negative weights, or a bad heuristic?"
                )
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                # heuristic = GEODESIC.inverse((lon, lat), dest_lonlat)[0, 0] / 3
                push(fringe, (vu_dist + heuristic, vu_dist, next(c), u))
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
            [
                TABLES["nodes"].c.node_id,
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
            ]
        )
        .where(
            (TABLES["nodes"].c.longitude != None) & (TABLES["nodes"].c.latitude != None)
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


nashville_lonlat = (-86.9759, 36.0346)
natchez_lonlat = (-91.36892, 31.54543)
# natchez_lonlat = (-89.36892, 33.31543)
nashville = find_node(*nashville_lonlat)
natchez = find_node(*natchez_lonlat)
print(nashville, natchez)

# Plotting
proj = ccrs.PlateCarree()
ax = plt.axes(projection=proj)

print("Plotting parkway…")
natchez_trace_parkway = shapely.wkb.load(open("./natchez_trace.wkb", "rb"))
ax.add_geometries(
    natchez_trace_parkway,
    ccrs.PlateCarree(),
    facecolor="green",
    alpha=0.8,
    edgecolor=[0, 0.5, 0, 1],
    zorder=32,
    linewidth=4,
)

print("Nashville to Natchez…")
pred = {nashville: (None, "self-loop")}
d = distances_from_focus(
    nashville, natchez, pred=pred, filter_nodes=True, filter_sources={"grid"}
)
print(d[natchez])
backtrack = natchez
lons = [natchez_lonlat[0]]
lats = [natchez_lonlat[1]]
sources = []
print("Backtrack…")
while pred.get(backtrack):
    backtrack, source = pred[backtrack]
    lonlat = DATABASE.execute(
        select(
            [
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
            ]
        ).where(TABLES["nodes"].c.node_id == backtrack)
    ).fetchone()
    if lonlat is None:
        break
    lons.append(lonlat[0])
    lats.append(lonlat[1])
    sources.append(source)

bbox = (min(lons), max(lons), min(lats), max(lats))
ax.plot(lons, lats, zorder=30, linewidth=4, c="blue")

print("Natchez to Nashville…")
rpred = {natchez: (None, "self-loop")}
rd = distances_from_focus(natchez, nashville, pred=rpred)
print(rd[nashville])
rbacktrack = nashville
rlons = [nashville_lonlat[0]]
rlats = [nashville_lonlat[1]]
rvisited = []
rsources = []
print("Backtrack…")
while rpred.get(rbacktrack):
    rbacktrack, source = rpred[rbacktrack]
    rvisited.append(rbacktrack)
    lonlat = DATABASE.execute(
        select(
            [
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
            ]
        ).where(TABLES["nodes"].c.node_id == rbacktrack)
    ).fetchone()
    if lonlat is None:
        break
    if len(rsources) >= 100:
        print(rbacktrack)
        breakpoint()
        break
    rlons.append(lonlat[0])
    rlats.append(lonlat[1])
    rsources.append(source)


print("Set up plot…")
bbox = (
    min(bbox[0], *rlons),
    max(bbox[1], *rlons),
    min(bbox[2], *rlats),
    max(bbox[3], *rlats),
)
ax.plot(rlons, rlats, zorder=31, linewidth=3, c="red")

ax.set_extent(bbox)  # x0, x1, y0, y1

random = ListedColormap(
    [(1, 1, 1, 0)] + [numpy.random.random(3) * 0.3 + 0.7 for _ in range(2 ** 16 - 1)],
    name="random",
)

print("Plot voronoi…")
for tile in tqdm(itertools.product(["N"], [10, 30], ["W"], [90, 120])):
    fname_v = "voronoi-{:s}{:d}{:s}{:d}.tif".format(*tile)
    fname_d = "distances-{:s}{:d}{:s}{:d}.tif".format(*tile)
    try:
        data = rasterio.open(fname_v, "r").read(1)[500:-500, 500:-500]
    except rasterio.errors.RasterioIOError:
        continue
    try:
        hs = 0
        for i in range(1, 9):
            hs = rasterio.open(fname_d, "r").read(i)[500:-500, 500:-500] + hs
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
    h = ax.imshow(
        numpy.minimum(hs, 7200),
        transform=proj,
        extent=(west, west + 30, south, south + 20),
        interpolation="nearest",
        cmap=cm.gist_yarg,
        alpha=0.4,
        vmax=7200,
    )
    h.set_zorder(2)

print("Background of plot…")
# shape_feature = ShapelyFeature(LAND, ccrs.PlateCarree(), facecolor="none", edgecolor="blue", lw=1)
ax.add_geometries(
    shpreader.Reader(
        shpreader.natural_earth(resolution="10m", category="physical", name="lakes")
    ).geometries(),
    ccrs.PlateCarree(),
    facecolor="cyan",
    alpha=0.8,
    edgecolor="none",
    zorder=4,
)

ax.add_geometries(
    shpreader.Reader(
        shpreader.natural_earth(resolution="10m", category="physical", name="ocean")
    ).geometries(),
    ccrs.PlateCarree(),
    facecolor="cyan",
    edgecolor="none",
    zorder=4,
)

hexagon_id, x, y, hx, hy, coast = zip(
    *DATABASE.execute(
        select(
            [
                TABLES["nodes"].c.node_id,
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
                TABLES["nodes"].c.h3longitude,
                TABLES["nodes"].c.h3latitude,
                TABLES["nodes"].c.coastal,
            ]
        ).where(TABLES["nodes"].c.node_id > 100000000)
    ).fetchall()
)

encountered = [(h in d, h in rd) for h in hexagon_id]


print("Nodes…")
ax.scatter(
    hx,
    hy,
    marker=".",
    edgecolors="none",
    facecolors=[
        {
            (False, False): [0, 0, 0, 1],
            (False, True): [0, 0, 1, 1],
            (True, True): [1, 0, 1, 1],
            (True, False): [1, 0, 0, 1],
        }[c]
        for c in encountered
    ],
    s=7,
    zorder=3,
)
ax.scatter(
    x,
    y,
    marker="*",
    facecolors="k",
    s=18,
    zorder=33,
)
node_id, x, y, coast = zip(
    *DATABASE.execute(
        select(
            [
                TABLES["nodes"].c.node_id,
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
                TABLES["nodes"].c.coastal,
            ]
        ).where(
            (TABLES["nodes"].c.node_id <= 100000000)
            & (TABLES["nodes"].c.node_id.in_(rvisited))
        )
    ).fetchall(),
)
ax.scatter(
    x,
    y,
    c="teal",
    marker="+",
    s=25,
    zorder=34,
)
ax.gridlines(draw_labels=True)


def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return int((lon + 180) / 6) + 1


from matplotlib import patheffects


def scale_bar(
    ax, proj, length, location=(0.5, 0.05), linewidth=3, units="km", m_per_unit=1000
):
    """

    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    # Projection in metres
    utm = ccrs.UTM(utm_from_lon((x0 + x1) / 2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit / 2, sbcx + length * m_per_unit / 2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(
        bar_xs,
        [sbcy, sbcy],
        transform=utm,
        color="k",
        linewidth=linewidth,
        path_effects=buffer,
    )
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(
        sbcx,
        sbcy,
        str(length) + " " + units,
        transform=utm,
        horizontalalignment="center",
        verticalalignment="bottom",
        path_effects=buffer,
        zorder=200,
    )
    right = x0 + (x1 - x0) * 0.95
    # Plot the N arrow
    t1 = ax.text(
        right,
        sbcy,
        "\u25B2\nN",
        transform=utm,
        horizontalalignment="center",
        verticalalignment="bottom",
        path_effects=buffer,
        zorder=200,
    )
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(
        bar_xs, [sbcy, sbcy], transform=utm, color="k", linewidth=linewidth, zorder=3
    )


print("Add scale bar…")
scale_bar(ax, ccrs.Mercator(), 100)  # 100 km scale bar

print("Show…")
plt.show()
