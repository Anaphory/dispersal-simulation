from itertools import count
from heapq import heappush as push, heappop as pop

import numpy
from tqdm import tqdm
from sqlalchemy import select
from sqlalchemy.dialects.sqlite import insert

import shapely
import shapely.geometry as sgeom
from shapely.prepared import prep
from h3.api import basic_int as h3
import cartopy.geodesic as geodesic
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader

from raster_data import *
from database import db
from ecoregions import TC

from matplotlib import pyplot as plt

# Define some constants
GEODESIC: geodesic.Geodesic = geodesic.Geodesic()

BBOX = sgeom.box(-168.956, -55.852, -17.556, 83.696)  # minx, miny, maxx, maxy

ALL_LAND = unary_union(
    [
        record.geometry
        for record in shpreader.Reader(
            shpreader.natural_earth(resolution="10m", category="physical", name="land")
        ).records()
        if record.attributes.get("featurecla") != "Null island"
    ]
).difference(
    unary_union(
        list(
            shpreader.Reader(
                shpreader.natural_earth(
                    resolution="10m", category="physical", name="lakes"
                )
            ).geometries()
        )
    )
)
LAND = BBOX.intersection(ALL_LAND)
PLAND = prep(LAND)

DATABASE, TABLES = db()


# =================================
# General geometry helper functions
# =================================
def as_point(h3_index):
    y, x = h3.h3_to_geo(index)
    return sgeom.Point(x, y)


# =========================
# Travelling time functions
# =========================
def navigation_speed(slope: float) -> float:
    """Using the slope in %, calculate the navigation speed in m/s

    This function calculates the off-road navigation speed (for male cadets in
    forested areas with navigational targets every 20 minutes) following
    [@irmischer2018measuring]. Like their formula, slope is in % and speed in
    m/s.

    > [T]he fastest off-road navigation speed was 0.78 m/s for males […] with a
    > peak at −2%.

    >>> navigation_speed(-2.)
    0.78
    >>> navigation_speed(-2.) > navigation_speed(-1.5)
    True
    >>> navigation_speed(-2.) > navigation_speed(-2.5)
    True

    """
    return 0.11 + 0.67 * numpy.exp(-((slope + 2.0) ** 2) / 1800.0)


def find_coast_hexagons():
    coastal = set()
    # Coast resolution is 10m and coasts are not known to be straight, so we
    # can expect that every coastal hexagon contains at least one of the
    # coastline polygon coordinates.
    for geom in LAND.boundary.geoms:
        for x, y in geom.coords:
            coastal.add(h3.geo_to_h3(y, x, 5))
    return coastal


# =================
# Raster operations
# =================
def prep_tile(
    lon: float,
    lat: float,
):
    print("Preparing rasters around ({:}, {:}):".format(lon, lat))
    elevation_file = gmted_tile_from_geocoordinates(lon, lat)
    m_trafo = elevation_file.transform
    height, width = elevation_file.shape
    elevation = numpy.full((height + 1000, width + 1000), -100, int)
    ecoregions = numpy.full((height + 1000, width + 1000), 999, int)
    elevation[500:-500, 500:-500] = elevation_file.read(1)
    ecoregions[500:-500, 500:-500] = ecoregion_tile_from_geocoordinates(lon, lat).read(
        1
    )
    print("Loading adjacent data…")
    try:
        elevation[:500, :500] = (
            gmted_tile_from_geocoordinates(lon - 30, lat + 20)
        ).read(1)[-500:, -500:]
        ecoregions[:500, :500] = ecoregion_tile_from_geocoordinates(
            lon - 30, lat + 20
        ).read(1)[-500:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[:500, 500:-500] = (
            gmted_tile_from_geocoordinates(lon, lat + 20)
        ).read(1)[-500:, :]
        ecoregions[:500, 500:-500] = ecoregion_tile_from_geocoordinates(
            lon, lat + 20
        ).read(1)[-500:, :]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[:500, -500:] = (
            gmted_tile_from_geocoordinates(lon + 30, lat + 20)
        ).read(1)[-500:, :500]
        ecoregions[:500, -500:] = ecoregion_tile_from_geocoordinates(
            lon + 30, lat + 20
        ).read(1)[-500:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[500:-500, :500] = (
            gmted_tile_from_geocoordinates(lon - 30, lat)
        ).read(1)[:, -500:]
        ecoregions[500:-500, :500] = ecoregion_tile_from_geocoordinates(
            lon - 30, lat
        ).read(1)[:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[500:-500, -500:] = (
            gmted_tile_from_geocoordinates(lon + 30, lat)
        ).read(1)[:, :500]
        ecoregions[500:-500, -500:] = ecoregion_tile_from_geocoordinates(
            lon + 30, lat
        ).read(1)[:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, :500] = (
            gmted_tile_from_geocoordinates(lon - 30, lat - 20)
        ).read(1)[:500, -500:]
        ecoregions[-500:, :500] = ecoregion_tile_from_geocoordinates(
            lon - 30, lat - 20
        ).read(1)[:500, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, 500:-500] = (
            gmted_tile_from_geocoordinates(lon, lat - 20)
        ).read(1)[:500, :]
        ecoregions[-500:, 500:-500] = ecoregion_tile_from_geocoordinates(
            lon, lat - 20
        ).read(1)[:500, :]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, -500:] = (
            gmted_tile_from_geocoordinates(lon + 30, lat - 20)
        ).read(1)[:500, :500]
        ecoregions[-500:, -500:] = ecoregion_tile_from_geocoordinates(
            lon + 30, lat - 20
        ).read(1)[:500, :500]
    except rasterio.RasterioIOError:
        pass

    transform = rasterio.Affine(
        m_trafo.a,
        0,
        m_trafo.c - 500 * m_trafo.a,
        0,
        m_trafo.e,
        m_trafo.f - 500 * m_trafo.e,
    )
    return elevation, ecoregions, transform


def all_pairwise_distances(
    elevation: numpy.array,
    transform: rasterio.Affine,
    terrain_coefficients: numpy.array,
) -> t.Dict[t.Tuple[int, int], numpy.array]:
    d_n, d_e, d_ne = [], [], []
    for y in range(1, len(elevation) + 1):
        (lon0, lat0) = transform * (0, y)
        (lon1, lat1) = transform * (1, y - 1)

        d = GEODESIC.inverse((lon0, lat0), [(lon0, lat1), (lon1, lat0), (lon1, lat1)])
        d_n.append(d[0, 0])
        d_e.append(d[1, 0])
        d_ne.append(d[2, 0])
    distance_to_north = numpy.array(d_n)[:-1]
    slope_to_north = (
        100 * (elevation[1:, :] - elevation[:-1, :]) / distance_to_north[:, None]
    )
    tc_to_north = (terrain_coefficients[1:, :] + terrain_coefficients[:-1, :]) / 2
    north = distance_to_north[:, None] / (
        navigation_speed(slope_to_north) * tc_to_north
    )
    south = distance_to_north[:, None] / (
        navigation_speed(-slope_to_north) * tc_to_north
    )
    del distance_to_north, slope_to_north, tc_to_north

    distance_to_east = numpy.array(d_e)
    slope_to_east = (
        100 * (elevation[:, 1:] - elevation[:, :-1]) / distance_to_east[:, None]
    )
    tc_to_east = (terrain_coefficients[:, 1:] + terrain_coefficients[:, :-1]) / 2
    east = distance_to_east[:, None] / (navigation_speed(slope_to_east) * tc_to_east)
    west = distance_to_east[:, None] / (navigation_speed(-slope_to_east) * tc_to_east)
    del distance_to_east, slope_to_east, tc_to_east

    distance_to_northeast = numpy.array(d_ne)[:-1]
    slope_to_northeast = (
        100 * (elevation[1:, 1:] - elevation[:-1, :-1]) / distance_to_northeast[:, None]
    )
    tc_to_northeast = (
        terrain_coefficients[1:, 1:] + terrain_coefficients[:-1, :-1]
    ) / 2
    northeast = distance_to_northeast[:, None] / (
        navigation_speed(slope_to_northeast) * tc_to_northeast
    )
    southwest = distance_to_northeast[:, None] / (
        navigation_speed(-slope_to_northeast) * tc_to_northeast
    )
    del distance_to_northeast, slope_to_northeast, tc_to_northeast
    distance_to_northwest = numpy.array(d_ne)[:-1]
    slope_to_northwest = (
        100 * (elevation[1:, :-1] - elevation[:-1, 1:]) / distance_to_northwest[:, None]
    )

    tc_to_northwest = (
        terrain_coefficients[1:, :-1] + terrain_coefficients[:-1, 1:]
    ) / 2

    southeast = distance_to_northwest[:, None] / (
        navigation_speed(-slope_to_northwest) * tc_to_northwest
    )
    northwest = distance_to_northwest[:, None] / (
        navigation_speed(slope_to_northwest) * tc_to_northwest
    )
    del distance_to_northwest, slope_to_northwest, tc_to_northwest

    return {
        (-1, 0): north,
        (-1, 1): northeast,
        (0, 1): east,
        (1, 1): southeast,
        (1, 0): south,
        (1, -1): southwest,
        (0, -1): west,
        (-1, -1): northwest,
    }


def distances_from_focus(
    source: t.Tuple[int, int],
    destinations: t.Optional[t.Set[t.Tuple[int, int]]],
    distance_by_direction: t.Dict[t.Tuple[int, int], numpy.array],
    pred: t.Optional,
) -> numpy.array:
    # Dijkstra's algorithm, adapted for our purposes from
    # networkx/algorithms/shortest_paths/weighted.html
    d = distance_by_direction[0, 1]
    dist: numpy.array = numpy.full((d.shape[0], d.shape[1] + 1), numpy.nan, dtype=float)
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
            r = min(r0, r1)
            c = min(c0, c1)
            if 0 <= r < d.shape[0] and 0 <= c < d.shape[1]:
                yield (r1, c1), d[r, c]

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


# ================
# Process hexagons
# ================
def find_land_hexagons():
    hexagons = find_coast_hexagons()
    for poly in LAND.geoms:
        d = shapely.geometry.mapping(poly)
        hexagons |= h3.polyfill_geojson(d, 5)
    return hexagons

plotting = False
if __name__ == "__main__":
    print("Plotting on")
    plotting = True

def core_point(hexbin, distance_by_direction, transform):
    def rowcol(latlon):
        lat, lon = latlon
        if lon > 170:
            # FIXME: We can and need to do this because we are working on the
            # Americas and the Americas only. The generic solution is more
            # difficult.
            lon = lon - 360
        col, row = ~transform * (lon, lat)
        return int(row), int(col)

    points = [rowcol(latlon) for latlon in h3.h3_to_geo_boundary(hexbin)]
    rmin = min(r for r, c in points) - 1
    rmax = max(r for r, c in points) + 1
    cmin = min(c for r, c in points) - 1
    cmax = max(c for r, c in points) + 1
    points = [(r - rmin, c - cmin) for r, c in points]

    dist = {
        (n, e): d[
            rmin - min(n, 0) : rmax - max(0, n), cmin - min(e, 0) : cmax - max(0, e)
        ]
        for (n, e), d in distance_by_direction.items()
    }

    border = []
    for r in range(rmin, rmax):
        for c in range(cmin, cmax):
            lon, lat = transform * (c, r)
            if h3.geo_to_h3(lat, lon, 5) == hexbin:
                x, y = transform * (c + 1, r)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((r - rmin, c - cmin))
                    continue
                x, y = transform * (c - 1, r)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((r - rmin, c - cmin))
                    continue
                x, y = transform * (c, r + 1)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((r - rmin, c - cmin))
                    continue
                x, y = transform * (c, r - 1)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((r - rmin, c - cmin))
                    continue

    if plotting:
        plt.imshow(dist[1, 0].T, extent=(0, cmax-cmin+1, rmax-rmin, 0))

    c = t.Counter()
    for r0, c0 in border:
        pred = {(r0, c0): None}
        all_dist = distances_from_focus((r0, c0), set(border), dist, pred=pred)
        for b1 in border:
            n = b1
            while pred[n]:
                if plotting:
                    plt.plot([n[1], pred[n][1]], [n[0], pred[n][0]], c="k")
                n = pred[n]
                c[n] += 1

    (r0, c0), centrality = c.most_common(1)[0]
    lon, lat = transform * (c0 + cmin, r0 + rmin)

    if plotting:
        rlat, rlon = h3.h3_to_geo(hexbin)
        print(f"Core node at ({lon}, {lat}). [Actual center at ({rlon}, {rlat}).]")
        rs, cs = zip(*border)
        plt.scatter(cs, rs)
        plt.scatter(c0, r0)
        rc, rr = ~transform * (rlon, rlat)
        plt.scatter(rc - cmin, rr - rmin)
        plt.show()

    return lon, lat


COAST = find_coast_hexagons()

hexes_by_tile = sorted(
    find_land_hexagons()
    - {n for n, in DATABASE.execute(select(TABLES["nodes"].c.node_id)).fetchall()},
    key=lambda hex: hash(tile_from_geocoordinates(*reversed(h3.h3_to_geo(hex)))),
)
tile = None
for hexagon in tqdm(hexes_by_tile):
    hlat, hlon = h3.h3_to_geo(hexagon)

    if tile != tile_from_geocoordinates(hlon, hlat):
        try:
            elevation, ecoregions, transform = prep_tile(hlon, hlat)
        except rasterio.errors.RasterioIOError:
            continue
        tile = tile_from_geocoordinates(hlon, hlat)
        terrain_coefficient_raster = TC[ecoregions]
        distance_by_direction = all_pairwise_distances(
            elevation, transform, terrain_coefficient_raster
        )

    lon, lat = core_point(hexagon, distance_by_direction, transform)
    DATABASE.execute(
        insert(TABLES["nodes"])
        .values(
            node_id=hexagon,
            longitude=lon,
            latitude=lat,
            h3longitude=hlon,
            h3latitude=hlat,
            coastal=(hexagon in COAST),
        )
        .on_conflict_do_nothing()
    )


if __name__ == "__main__":
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
    x, y, coast = zip(
        *DATABASE.execute(
            select(
                TABLES["nodes"].c.longitude,
                TABLES["nodes"].c.latitude,
                TABLES["nodes"].c.coastal,
            )
        ).fetchall()
    )
    shape_feature = ShapelyFeature(
        LAND, ccrs.PlateCarree(), facecolor="none", edgecolor="blue", lw=1
    )
    ax.scatter(x, y, c=coast)
    ax.add_feature(shape_feature)
    plt.show()
