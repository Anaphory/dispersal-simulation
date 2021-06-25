from itertools import count
import typing as t
from heapq import heappush as push, heappop as pop

from sqlalchemy import func
import numpy
from numpy import pi, cos, inf
from sqlalchemy import select

import cartopy.geodesic as geodesic

from database import db
from earth import LAND, BBOX, GEODESIC

DATABASE, TABLES = db()


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
            continue  # already searched this node.
        dist[spot] = d
        if spot == destination:
            break

        for u, cost in DATABASE.execute(
            select(
                TABLES["edges"].c.node2,
                TABLES["edges"].c.travel_time,
            ).where(TABLES["edges"].c.node1 == spot, *filter)
        ):
            vu_dist = dist[spot] + cost
            if numpy.isfinite(dist[u]):
                if vu_dist < dist[u]:
                    raise ValueError(
                        "Contradictory paths found. Do you have negative weights?"
                    )
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if pred is not None:
                    pred[u] = spot
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
pred = {nashville: None}
distances_from_focus(nashville, natchez, pred=pred)
backtrack = natchez
while pred.get(backtrack):
    print(backtrack)
    backtrack = pred(backtrack)

import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.feature import ShapelyFeature
from matplotlib import pyplot as plt

proj = ccrs.PlateCarree()
ax = plt.axes(projection=proj)
ax.set_extent(
    (BBOX.bounds[0], BBOX.bounds[2], BBOX.bounds[1], BBOX.bounds[3])
)  # x0, x1, y0, y1
ax.stock_img()

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
c = numpy.array(coast)
ax.scatter(hx, hy, c=c, marker="*")
ax.scatter(x, y, c=c, marker=".")

node_id, x, y, coast = zip(
    *DATABASE.execute(
        select(
            TABLES["nodes"].c.node_id,
            TABLES["nodes"].c.longitude,
            TABLES["nodes"].c.latitude,
            TABLES["nodes"].c.coastal,
        ).where(TABLES["nodes"].c.node_id <= 100000000)
    ).fetchall()
)
c = numpy.array(coast)
ax.scatter(x, y, c=c, marker="+")


shape_feature = ShapelyFeature(
    LAND, ccrs.PlateCarree(), facecolor="none", edgecolor="blue", lw=1
)
ax.add_feature(shape_feature)
plt.show()
