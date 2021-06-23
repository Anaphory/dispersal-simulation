import shapely
import shapely.geometry as sgeom
from shapely.prepared import prep
from h3.api import basic_int as h3
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader


from raster_data import *
from database import db

BBOX = sgeom.box(-168.956, -55.852, -17.556, 83.696) # minx, miny, maxx, maxy

all_land = unary_union([
    record.geometry
    for record in shpreader.Reader(shpreader.natural_earth(
    resolution='10m', category='physical', name='land')).records()
    if record.attributes.get('featurecla') != "Null island"])

our_land = BBOX.intersection(all_land)

PLAND = prep(our_land)

def main():
    database, tables = db()

def find_coast_hexagons():
    coastal = set()
    # Coast resolution is 10m and coasts are not known to be straight, so we
    # can expect that every coastal hexagon contains at least one of the
    # coastline polygon coordinates.
    for geom in our_land.boundary.geoms:
        for x, y in geom.coords:
            coastal.add(h3.geo_to_h3(y, x, 5))
    return coastal

def as_point(index):
    y, x = h3.h3_to_geo(index)
    return sgeom.Point(x, y)

def prep_tile(
        lon: float, lat: float,
):
    """Compute pairwise distances and ecoregion composition

    Given a digital elevation model as raster map and a matching ecoregions
    raster map, compute the pairwise distance to its 1-hex and 2-hex neighbors
    for every H3 address hex at standard resolution, as well as the approximate
    cover of that cell in terms of ecoregions, for all cells where that is
    possible.

    Distances are computed in gross hours of travel while navigating off-track,
    following [@irmischer2018measuring].

    Returns
    =======
    d: A mapping. d[h1][h2] is the distance, in hours, from the center of h1 to
        the center of h2.
    e: A mapping. d[h1][b] is the proportion of hex h1 covered by ecoregion b.

    """
    print("Working on hex around ({:}, {:}):".format(lon, lat))
    elevation_file = gmted_tile_from_geocoordinates(lon, lat)
    m_trafo = elevation_file.transform
    height, width = elevation_file.shape
    elevation = numpy.full((height + 1000, width + 1000), -100, int)
    ecoregions = numpy.full((height + 1000, width + 1000), 999, int)
    elevation[500:-500, 500:-500] = elevation_file.read(1)
    ecoregions[500:-500, 500:-500] = ecoregion_tile_from_geocoordinates(lon, lat).read(1)
    print("Loading adjacent data…")
    try:
        elevation[:500, :500] = (gmted_tile_from_geocoordinates(lon - 30, lat + 20)).read(1)[-500:, -500:]
        ecoregions[:500, :500] = ecoregion_tile_from_geocoordinates(lon - 30, lat + 20).read(1)[-500:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[:500, 500:-500] = (gmted_tile_from_geocoordinates(lon, lat + 20)).read(1)[-500:, :]
        ecoregions[:500, 500:-500] = ecoregion_tile_from_geocoordinates(lon, lat + 20).read(1)[-500:, :]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[:500, -500:] = (gmted_tile_from_geocoordinates(lon + 30, lat + 20)).read(1)[-500:, :500]
        ecoregions[:500, -500:] = ecoregion_tile_from_geocoordinates(lon + 30, lat + 20).read(1)[-500:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[500:-500, :500] = (gmted_tile_from_geocoordinates(lon - 30, lat)).read(1)[:, -500:]
        ecoregions[500:-500, :500] = ecoregion_tile_from_geocoordinates(lon - 30, lat).read(1)[:, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[500:-500, -500:] = (gmted_tile_from_geocoordinates(lon + 30, lat)).read(1)[:, :500]
        ecoregions[500:-500, -500:] = ecoregion_tile_from_geocoordinates(lon + 30, lat).read(1)[:, :500]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, :500] = (gmted_tile_from_geocoordinates(lon - 30, lat - 20)).read(1)[:500, -500:]
        ecoregions[-500:, :500] = ecoregion_tile_from_geocoordinates(lon - 30, lat - 20).read(1)[:500, -500:]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, 500:-500] = (gmted_tile_from_geocoordinates(lon, lat - 20)).read(1)[:500, :]
        ecoregions[-500:, 500:-500] = ecoregion_tile_from_geocoordinates(lon, lat - 20).read(1)[:500, :]
    except rasterio.RasterioIOError:
        pass
    try:
        elevation[-500:, -500:] = (gmted_tile_from_geocoordinates(lon + 30, lat - 20)).read(1)[:500, :500]
        ecoregions[-500:, -500:] = ecoregion_tile_from_geocoordinates(lon + 30, lat - 20).read(1)[:500, :500]
    except rasterio.RasterioIOError:
        pass

    print("Computing hex extents…")
    transform = rasterio.Affine(m_trafo.a, 0, m_trafo.c - 500 * m_trafo.a,
                                0, m_trafo.e, m_trafo.f - 500 * m_trafo.e)
    return elevation, ecoregions, transform

def find_land_hexagons():
    hexagons = find_coast_hexagons()
    for poly in our_land.geoms:
        d = shapely.geometry.mapping(poly)
        hexagons |= h3.polyfill_geojson(d, 5)
    return hexagons

from ecoregions import TC

def core_point(hexbin):
    lat, lon = h3.h3_to_geo(hexbin)
    elevation, ecoregions, transform = prep_tile(lon, lat)
    terrain_coefficient_raster = TC[ecoregions]
    distance_by_direction = all_pairwise_distances(
        elevation, transform, terrain_coefficient_raster)

    def rowcol(latlon):
        lat, lon = latlon
        if lon > 170:
            # FIXME: We can and need to do this because we are working on the
            # Americas and the Americas only. The generic solution is more
            # difficult.
            lon = lon - 360
        col, row = ~transform * (lon, lat)
        return int(row), int(col)

    points = [rowcol(latlon)
              for latlon in h3.h3_to_geo_boundary(hexbin)]
    rmin = min(r for r, c in points)
    rmax = max(r for r, c in points) + 1
    cmin = min(c for r, c in points)
    cmax = max(c for r, c in points) + 1

    dist = {(n, e): d[rmin-min(n, 0):rmax  - max(0, n),
                      cmin-min(e, 0):cmax  - max(0, e)]
            for (n, e), d in distance_by_direction.items()}

    border = []
    for r in range(rmin, rmax):
        for c in range(cmin, cmax):
            lon, lat = transform * (c, r)
            if h3.geo_to_h3(lat, lon, 5) == hexbin:
                x, y = transform * (c + 1, r)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((c, r))
                    continue
                x, y = transform * (c - 1, r)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((c, r))
                    continue
                x, y = transform * (c, r + 1)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((c, r))
                    continue
                x, y = transform * (c, r - 1)
                if h3.geo_to_h3(y, x, 5) != hexbin:
                    border.append((c, r))
                    continue

    c = t.Counter()
    for r0, c0 in border:
        pred = {(r0, c0): None}
        all_dist = distances_from_focus((r0, c0), set(border), dist, pred=pred)
        for b1 in border:
            n = b1
            while pred[n]:
                n = pred[n]
                c[n] += 1
    (r0, c0), centrality = c.most_common(1)[0]
    center[hexbin] = (r0 + rmin, c0 + cmin)
    lon, lat = transform * (c0 + cmin, r0 + rmin)
    rlat, rlon = h3.h3_to_geo(hexbin)
    print(f"Centalic node at ({lon}, {lat}). [Actual center at ({rlon}, {rlat}).]")
    return (c0 + cmin, r0 + rmin), transform



if __name__=="__main__":
    main()

    import cartopy.crs as ccrs
    import cartopy.feature as cf
    from cartopy.feature import ShapelyFeature
    from matplotlib import pyplot as plt
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.set_extent((BBOX.bounds[0], BBOX.bounds[2], BBOX.bounds[1], BBOX.bounds[3])) # x0, x1, y0, y1
    ax.stock_img()
    shape_feature = ShapelyFeature(our_land, ccrs.PlateCarree(), facecolor='none', edgecolor='blue', lw=1)
    y, x = zip(*[h3.h3_to_geo(h)    for h in     find_land_hexagons()])
    ax.scatter(x, y)
    ax.add_feature(shape_feature)
    plt.show()
