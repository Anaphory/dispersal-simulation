#!/home/cluster/gkaipi/.pyenv/shims/python
import json
import itertools
import typing as t
import sys
from itertools import count
from heapq import heappush as push, heappop as pop

import numpy
from more_itertools import windowed
from tqdm import tqdm
from sqlalchemy import func
from sqlalchemy import select
from sqlalchemy.dialects.sqlite import insert

import shapely
import shapely.geometry as sgeom
from h3.api import basic_int as h3

import rasterio
from raster_data import ecoregion_tile, boundingbox_from_tile, Tile, RowCol
from database import db
from earth import LAND, GEODESIC, LonLat, DEFINITELY_INLAND
from by_river import process_rivers
from by_sea import distance_by_sea
from stitch import merge_tile

DATABASE, TABLES = db()
RESOLUTION = 5


# ============
# Type aliases
# ============
H3Index = int


# =================================
# General geometry helper functions
# =================================
def as_point(h3_index: H3Index) -> sgeom.Point:
    y, x = h3.h3_to_geo(h3_index)
    return sgeom.Point(x, y)


def find_coast_hexagons() -> t.Set[H3Index]:
    coastal = set()
    # Coast resolution is 10m and coasts are not known to be straight, so we
    # can expect that every coastal hexagon contains at least one of the
    # coastline polygon coordinates.
    for geom in LAND.boundary.geoms:
        for x, y in geom.coords:
            coastal.add(h3.geo_to_h3(y, x, 5))
    return coastal


try:
    COAST = set(json.load(open("COAST.json")))
except (FileNotFoundError, json.JSONDecodeError):
    COAST = find_coast_hexagons()
    json.dump(list(COAST), open("COAST.json", "w"))


# ==========================
# Distance raster operations
# ==========================


def load_distances(tile: Tile):
    fname = "distances-{:s}{:d}{:s}{:d}.tif".format(*tile)

    distance_raster = rasterio.open(fname)
    return {
        (-1, 0): distance_raster.read(1),
        (-1, 1): distance_raster.read(2),
        (0, 1): distance_raster.read(3),
        (1, 1): distance_raster.read(4),
        (1, 0): distance_raster.read(5),
        (1, -1): distance_raster.read(6),
        (0, -1): distance_raster.read(7),
        (-1, -1): distance_raster.read(8),
    }, distance_raster.transform


def distances_from_focus(
    source: t.Tuple[int, int],
    destinations: t.Optional[t.Set[RowCol]],
    distance_by_direction: t.Dict[RowCol, numpy.ndarray],
    pred: t.Optional[t.Dict[RowCol, RowCol]] = None,
) -> numpy.ndarray:
    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    d = distance_by_direction[0, 1]
    dist: numpy.array = numpy.full(d.shape, numpy.inf, dtype=float)
    seen: t.Dict[t.Tuple[int, int], float] = {source: 0.0}
    c = count()
    # use heapq with (distance, label) tuples
    fringe: t.List[t.Tuple[float, int, t.Tuple[int, int]]] = []
    push(fringe, (0, next(c), source))

    def moore_neighbors(
        r0: int, c0: int
    ) -> t.Iterable[t.Tuple[t.Tuple[int, int], float]]:
        for (r, c), d in distance_by_direction.items():
            r1, c1 = r0 + r, c0 + c
            if 0 <= r1 < d.shape[0] and 0 <= c1 < d.shape[1]:
                yield (r1, c1), d[r0, c0]

    for u, cost in moore_neighbors(*source):
        if not numpy.isfinite(cost):
            r0, c0 = source[0] - u[0], source[1] - u[1]
            # We don't have the actual distance away from the river, so we
            # assume the terrain is somewhat homogenous around it.
            cost = max(
                distance_by_direction[r0, c0][u], distance_by_direction[-r0, -c0][u]
            )
            if numpy.isfinite(cost):
                push(fringe, (cost, next(c), u))

    while fringe:
        (d, _, spot) = pop(fringe)
        if numpy.isfinite(dist[spot]):
            continue  # already searched this node.
        dist[spot] = d
        if destinations is not None and spot in destinations:
            destinations.remove(spot)
            if not destinations:
                break

        for u, cost in moore_neighbors(*spot):
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


def distances_trimmed(source, points, distance_by_direction, transform):
    rmin = min(r for r, c in points) - 2
    rmax = max(r for r, c in points) + 3
    cmin = min(c for r, c in points) - 2
    cmax = max(c for r, c in points) + 3
    points = [(r - rmin, c - cmin) for r, c in points]

    r0, c0 = source[0] - rmin, source[1] - cmin

    dist = {
        (n, e): d[rmin:rmax, cmin:cmax] for (n, e), d in distance_by_direction.items()
    }

    return rmin, cmin, distances_from_focus((r0, c0), points, dist, pred=None)


# ================
# Process hexagons
# ================
def find_land_hexagons():
    hexagons = set()
    for poly in LAND.geoms:
        d = shapely.geometry.mapping(poly)
        hexagons |= h3.polyfill_geojson(d, 5)
    return hexagons | COAST


try:
    ALL = set(json.load(open("ALL.json")))
except (FileNotFoundError, json.JSONDecodeError):
    ALL = find_land_hexagons()
    json.dump(list(ALL), open("ALL.json", "w"))


LONG_TIME = 12 * 3600


def core_point(hexbin, distance_by_direction, is_river, transform):
    def rowcol(latlon):
        lat, lon = latlon
        if lon > 170:
            # FIXME: We can and need to do this because we are working on the
            # Americas and the Americas only. The generic solution is more
            # difficult.
            lon = lon - 360
        col, row = ~transform * (lon, lat)
        return round(row), round(col)

    points = [rowcol(latlon) for latlon in h3.h3_to_geo_boundary(hexbin)]
    rmin = min(r for r, c in points) - 1
    rmax = max(r for r, c in points) + 2
    cmin = min(c for r, c in points) - 1
    cmax = max(c for r, c in points) + 2
    points = [(r - rmin, c - cmin) for r, c in points]

    dist = {
        # For the purposes of finding central locations, we cannot have so so
        # many locations each with the same infinite centrality, so assume a
        # maximum pixel distance of 8 hours.
        (n, e): numpy.minimum(d[rmin:rmax, cmin:cmax], LONG_TIME)
        for (n, e), d in distance_by_direction.items()
    }
    dist[0, -1] = dist[0, 1] = dist[0, 1] + dist[0, -1]
    dist[-1, -1] = dist[1, 1] = dist[1, 1] + dist[-1, -1]
    dist[1, -1] = dist[-1, 1] = dist[-1, 1] + dist[1, -1]
    dist[-1, 0] = dist[1, 0] = dist[1, 0] + dist[-1, 0]

    border = points[:]
    for i in range(-1, len(points) - 1):
        r0, c0 = points[i]
        r1, c1 = points[i + 1]
        border.append(
            (
                round(0.5 * (r0 + r1)),
                round(0.5 * (c0 + c1)),
            )
        )

    all_dist = 0
    for r0, c0 in border:
        all_dist = all_dist + distances_from_focus((r0, c0), None, dist)

    all_dist[is_river[rmin:rmax, cmin:cmax]] = numpy.inf
    (r0, c0) = numpy.unravel_index(numpy.argmin(all_dist), all_dist.shape)
    return (r0 + rmin, c0 + cmin)


def tile_core_points(tile: Tile, skip_existing=True, n=0):
    west, south, east, north = boundingbox_from_tile(tile)

    def tile_contains(hexagon: H3Index):
        lat, lon = h3.h3_to_geo(hexagon)
        return (west < lon < east) and (south < lat < north)

    hexes_by_tile = {hex for hex in ALL if tile_contains(hex)}
    if skip_existing:
        hexes_by_tile -= {
            n
            for n, in DATABASE.execute(
                select(TABLES["nodes"].c.node_id).where(
                    TABLES["nodes"].c.h3longitude >= west,
                    TABLES["nodes"].c.h3latitude >= south,
                    TABLES["nodes"].c.h3longitude < east,
                    TABLES["nodes"].c.h3latitude < north,
                )
            )
        }

    distance_by_direction, transform = load_distances(tile)
    rivers = rasterio.open("rivers-{:s}{:d}{:s}{:d}.tif".format(*tile))
    is_river = rivers.read(1).astype(bool)

    try:
        values = []
        for hexagon in tqdm(hexes_by_tile):
            hlat, hlon = h3.h3_to_geo(hexagon)
            row, col = core_point(hexagon, distance_by_direction, is_river, transform)
            lon, lat = transform * (col, row)
            n += 1
            values.append(
                {
                    "node_id": hexagon,
                    "longitude": lon,
                    "latitude": lat,
                    "h3longitude": hlon,
                    "h3latitude": hlat,
                    "coastal": (hexagon in COAST),
                    "short": n,
                }
            )
    finally:
        if values:
            DATABASE.execute(
                insert(TABLES["nodes"]).values(values).on_conflict_do_nothing()
            )
    return n


def voronoi_and_neighbor_distances(tile, skip_existing=True):
    west, south, east, north = boundingbox_from_tile(tile)
    nodes = [
        (node, short, (lon, lat))
        for node, short, lon, lat in DATABASE.execute(
            select(
                TABLES["nodes"].c.node_id,
                TABLES["nodes"].c.short,
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
            ).where(
                TABLES["nodes"].c.longitude >= west,
                TABLES["nodes"].c.latitude >= south,
                TABLES["nodes"].c.longitude < east,
                TABLES["nodes"].c.latitude < north,
            )
        )
    ]

    distance_by_direction, transform = load_distances(tile)
    fname_v = "voronoi-{:s}{:d}{:s}{:d}.tif".format(*tile)
    fname_d = "min_distances-{:s}{:d}{:s}{:d}.tif".format(*tile)
    try:
        voronoi_allocation = rasterio.open(fname_v).read(1)
        min_distances = rasterio.open(fname_d).read(1)
    except rasterio.errors.RasterioIOError:
        height, width = distance_by_direction[1, 1].shape
        voronoi_allocation = numpy.zeros((height, width), dtype=numpy.uint32)
        min_distances = numpy.full(
            (height, width), numpy.inf, dtype=numpy.float32
        )

    def rowcol(lonlat: LonLat):
        nlon, nlat = lonlat
        if nlon > 170:
            # FIXME: We can and need to do this because we are working on the
            # Americas and the Americas only. The generic solution is more
            # difficult.
            nlon = nlon - 360
        ncol, nrow = ~transform * (nlon, nlat)
        return (round(nrow), round(ncol))

    try:
        for node, short, (lon, lat) in tqdm(nodes):
            if node > 100000000:
                # Make sure we have a small int to refer to this node!
                if short is None:
                    raise ValueError("Short ID of a hex node must not be None")
                # Node is an h3 index, not a river reach index
                environment = h3.k_ring(node, 2)
            else:
                # Node is a river reach end point
                environment = h3.k_ring(h3.geo_to_h3(lat, lon, 5), 1)
                environment.add(node)

            extremities = DATABASE.execute(
                select(
                    TABLES["nodes"].c.longitude,
                    TABLES["nodes"].c.latitude,
                    TABLES["nodes"].c.h3longitude,
                    TABLES["nodes"].c.h3latitude,
                ).where(
                    TABLES["nodes"].c.longitude != None,
                    TABLES["nodes"].c.latitude != None,
                    TABLES["nodes"].c.node_id.in_(environment),
                )
            ).fetchall()

            x, y, hx, hy = zip(*extremities)

            source = rowcol((lon, lat))
            nnodes = []
            nrowcol = []
            neighbors = DATABASE.execute(
                select(
                    TABLES["nodes"].c.node_id,
                    TABLES["nodes"].c.longitude,
                    TABLES["nodes"].c.latitude,
                ).where(
                    min(x) <= TABLES["nodes"].c.longitude,
                    TABLES["nodes"].c.longitude <= max(x),
                    min(y) <= TABLES["nodes"].c.latitude,
                    TABLES["nodes"].c.latitude <= max(y),
                )
            ).fetchall() + [
                # Add the hexagon neighbors as virtual targets, to have a
                # uniform breadth
                (None, hlon, hlat)
                for hlon, hlat in zip(hx, hy)
                if hlon is not None and hlat is not None
            ]

            if skip_existing:
                already_known = {
                    n
                    for n, in DATABASE.execute(
                        select(TABLES["edges"].c.node2).where(
                            TABLES["edges"].c.node1 == node,
                            TABLES["edges"].c.source == "grid",
                        )
                    )
                }
                to_be_known = {n for n, _, _ in neighbors} - {None}
                if already_known >= to_be_known and (
                    (node < 100000000)
                    or (
                        min_distances[source] == 0.0
                        and voronoi_allocation[source] == short
                    )
                ):
                    continue
                print(already_known - to_be_known, to_be_known - already_known)

            for nnode, nlon, nlat in neighbors:
                row, col = rowcol((nlon, nlat))
                if (
                    0 <= row < distance_by_direction[1, 1].shape[0]
                    and 0 <= col < distance_by_direction[1, 1].shape[1]
                ):
                    nrowcol.append((row, col))
                    nnodes.append(nnode)

            rmin, cmin, array = distances_trimmed(
                source, nrowcol, distance_by_direction, transform
            )

            values = [
                {
                    "node1": node,
                    "node2": nnode,
                    "source": "grid",
                    "travel_time": array[nr - rmin, nc - cmin],
                }
                for nnode, (nr, nc) in zip(nnodes, nrowcol)
                if nnode is not None
            ]
            DATABASE.execute(
                insert(TABLES["edges"]).values(values).on_conflict_do_nothing()
            )

            if node > 100000000:
                # Node is an h3 index, not a river reach index, so update the voronoi shapes around it
                rows, cols = array.shape
                voronoi_allocation[rmin : rmin + rows, cmin : cmin + cols][
                    min_distances[rmin : rmin + rows, cmin : cmin + cols] >= array
                ] = short
                min_distances[rmin : rmin + rows, cmin : cmin + cols] = numpy.fmin(
                    min_distances[rmin : rmin + rows, cmin : cmin + cols], array
                )
    finally:
        profile = rasterio.profiles.DefaultGTiffProfile()
        profile["height"] = voronoi_allocation.shape[0]
        profile["width"] = voronoi_allocation.shape[1]
        profile["transform"] = transform
        del profile["dtype"]
        profile["count"] = 1

        with rasterio.open(
            fname_v,
            "w",
            dtype=numpy.uint32,
            **profile,
        ) as dst:
            dst.write(voronoi_allocation, 1)
        with rasterio.open(
            fname_d,
            "w",
            dtype=numpy.float32,
            **profile,
        ) as dst:
            dst.write(min_distances, 1)


def measure_ecoregions(tile: Tile):
    fname_v = "voronoi-{:s}{:d}{:s}{:d}.tif".format(*tile)
    # Only read the non-overlapping part
    voronoi_allocation = rasterio.open(fname_v).read(
        1, window=((500, 5800 - 500), (500, 8200 - 500))
    )
    ecoregion_data = ecoregion_tile(tile)
    ecoregions = ecoregion_data.read(1)
    transform = ecoregion_data.transform
    areas = numpy.zeros(len(ecoregions))
    for y, _ in tqdm(enumerate(areas)):
        # Approximate grid cells by rectangles.
        (lon0, lat0) = transform * (0, y)
        (lon1, lat1) = transform * (1, y + 1)

        d = GEODESIC.inverse((lon0, lat0), [(lon0, lat1), (lon1, lat0)])
        areas[y] = d[0, 0] * d[0, 1]

    values = t.DefaultDict(t.Counter)
    for area, voronoi, eco in zip(areas, voronoi_allocation, ecoregions):
        for v, e in zip(voronoi, eco):
            if e != 999 and v != 0:
                values[int(v)][int(e)] += area
    return values


if __name__ == "__main__":
    if "reversed" in sys.argv and "from-mid" in sys.argv:

        def order(x):
            x = list(reversed(list(x)))
            m = len(x) // 2
            return x[m:] + x[:m]

    elif "reversed" in sys.argv:

        def order(x):
            return reversed(list(x))

    elif "from-mid" in sys.argv:

        def order(x):
            x = list(x)
            m = len(x) // 2
            return x[m:] + x[:m]

    else:

        def order(x):
            return list(x)

    if "rivers" in sys.argv:
        for tile in order(
            itertools.product(
                ["N", "S"],
                [10, 30, 50, 70],
                ["E", "W"],
                [0, 30, 60, 90, 120, 150, 180],
            )
        ):
            try:
                process_rivers(tile)
            except rasterio.errors.RasterioIOError:
                print("Tile {:s}{:d}{:s}{:d} not found.".format(*tile))

    if "core" in sys.argv:
        (n,) = DATABASE.execute(select(func.max(TABLES["nodes"].c.short))).fetchone()
        if n is None:
            n = 0
        else:
            n = (n // 140000 + 1) * 140000
        for tile in order(
            itertools.product(
                ["N", "S"],
                [10, 30, 50, 70],
                ["E", "W"],
                [0, 30, 60, 90, 120, 150, 180],
            )
        ):
            try:
                n = tile_core_points(tile, n=n)
            except rasterio.errors.RasterioIOError:
                print("Tile {:s}{:d}{:s}{:d} not found.".format(*tile))

    if "popdense" in sys.argv:
        # A modification of Tallavaara's R Markdown script. This can be run at
        # this point or later, it needs to be finished *and merged* before the
        # `populations` step.
        ...

    if "sea" in sys.argv:
        # Separate module. This can run in parallel to the following steps, it
        # only needs all core points to be defined, but does not interfere with
        # voronoi cells or grid distances.
        distance_by_sea(DEFINITELY_INLAND)

    if "voronoi" in sys.argv or "distances" in sys.argv:
        for tile in order(
            itertools.product(
                ["N", "S"], [10, 30, 50, 70], ["E", "W"], [0, 30, 60, 90, 120, 150, 180]
            )
        ):
            print("Working on tile {:s}{:d}{:s}{:d}…".format(*tile))
            try:
                voronoi_and_neighbor_distances(tile)
            except rasterio.errors.RasterioIOError:
                print("Tile {:s}{:d}{:s}{:d} not found.".format(*tile))

    if "stitch" in sys.argv:
        # Separate module, and quite fast. You could run it on every Voronoi
        # tile whose neighbors have been computed.
        for lon in order([-15, -45, -75, -105, -135, -165]):
            for lat in order([0, 20, 40, 60, 80, -20, -40, -60, -80]):
                try:
                    merge_tile(lon, lat)
                except rasterio.errors.RasterioIOError:
                    print("Stitching: tile for {:d}, {:d} not found.".format(lon, lat))

    if "areas" in sys.argv:
        # The other steps can be run per tile (with some artefacts on the tile
        # boundaries, if the neighboring tiles have not finished the previous
        # step), but this one really needs all voronoi computations to have
        # finished.
        values = t.DefaultDict(t.Counter)
        node_from_short = dict(
            list(
                DATABASE.execute(
                    select(
                        TABLES["nodes"].c.short,
                        TABLES["nodes"].c.node_id,
                    ).where(TABLES["nodes"].c.short != None)
                )
            )
        )
        for tile in order(
            itertools.product(
                ["N", "S"], [10, 30, 50, 70], ["E", "W"], [0, 30, 60, 90, 120, 150, 180]
            )
        ):
            print("Working on tile {:s}{:d}{:s}{:d}…".format(*tile))
            try:
                for short, counter in measure_ecoregions(tile).items():
                    values[node_from_short.get(short, -short)].update(counter)
            except rasterio.errors.RasterioIOError:
                print("Tile {:s}{:d}{:s}{:d} not found.".format(*tile))
            json.dump(values, open("areas.json", "w"), indent=2)
        all_data = [
            {"node": node, "ecoregion": ecoregion, "area": area / 1000000}
            for node, areas in values.items()
            for ecoregion, area in areas.items()
        ]
        for window in windowed(all_data, 300, 300):
            DATABASE.execute(insert(TABLES["ecology"]).values(window))

    if "populations" in sys.argv:
        # Needs `popdense` and `areas`.
        for node, ecoregion, area, popdensity in DATABASE.execute(
            select(
                TABLES["ecology"].c.node,
                TABLES["ecology"].c.ecoregion,
                TABLES["ecology"].c.area,
                TABLES["nodes"].c.popdensity,
            ).select_from(TABLES["ecology"].join(TABLES["nodes"]))
        ):
            DATABASE.execute(
                TABLES["ecology"]
                .update()
                .values(population_capacity=popdensity * area)
                .where(
                    (TABLES["ecology"].c.node == node)
                    & (TABLES["ecology"].c.ecoregion == ecoregion)
                )
            )
