import typing as t
from pathlib import Path

import more_itertools
import numpy

import rasterio
import rasterio.features
from h3.api import numpy_int as h3
import zipfile

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep

from ecoregions import ECOREGIONS, TC, RESOLUTION

import cartopy.geodesic as geodesic
GEODESIC: geodesic.Geodesic = geodesic.Geodesic()


def tile_from_geocoordinates(
        lon: float, lat: float
) -> t.Tuple[t.Literal["N", "S"],
             int,
             t.Literal["E", "W"],
             int]:
    """Turn a Longitude/Latitude (i.e. x, y) pair into a tile index

    The index describes the South West corner of the tile, which has a width of
    30° (anchored at the 0° meridian) and a height of 20° (centered at the
    equator, so anchored at ±10°). This kind of index is used eg. for the GMTED
    tiles.

    >>> tile_from_geocoordinates(15, 40)
    ('N', 30, 'E', 0)
    >>> tile_from_geocoordinates(-75, -20)
    ('S', 30, 'W', 90)
    >>> tile_from_geocoordinates(-179, -17)
    ('S', 30, 'W', 180)

    """
    southwest_corner_lon = int(lon // 30) * 30
    southwest_corner_lat = int((lat - 10) // 20) * 20 + 10
    ew: t.Literal["E", "W"] = "W" if southwest_corner_lon < 0 else "E"
    ns: t.Literal["N", "S"] = "S" if southwest_corner_lat < 0 else "N"
    return ns, abs(southwest_corner_lat), ew, abs(southwest_corner_lon)


def ecoregion_tile_from_geocoordinates(lon: float, lat: float) -> rasterio.DatasetReader:
    ns, lat0, ew, lon0 = tile_from_geocoordinates(lon, lat)
    ecoregions_path_t = "../ecoregions/ECOREGIONS-{0:02d}{1}{2:03d}{3}_20101117_gmted_med150.tif"
    return rasterio.open(ecoregions_path_t.format(lat0, ns.lower(), lon0, ew.lower()))


def gmted_tile_from_geocoordinates(lon: float, lat: float) -> rasterio.DatasetReader:
    """Path to the GMTED tile for given geocoordinates.

    The return type is explicitly a `str`, not a `pathlib.Path`, because it
    contains a GDAL virtual file system component.

    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010N30E000/
    >>> d1 = gmted_tile_from_geocoordinates(15, 40)
    >>> d1.name[-78:]
    'elevation/GMTED2010/GMTED2010N30E000_150.zip/30n000e_20101117_gmted_med150.tif'

    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010S30W090/
    >>> d2 = gmted_tile_from_geocoordinates(-75, -20)
    >>> d2.name[-78:]
    'elevation/GMTED2010/GMTED2010S30W090_150.zip/30s090w_20101117_gmted_med150.tif'

    """
    ns, lat, ew, lon = tile_from_geocoordinates(lon, lat)
    file = "{:02d}{:1s}{:03d}{:1s}".format(
        lat, ns.lower(), lon, ew.lower())
    tile = "GMTED2010{:1s}{:02d}{:1s}{:03d}".format(
        ns.upper(), lat, ew.upper(), lon)
    path = (Path(__file__).absolute().parent.parent /
            f"elevation/GMTED2010/{tile:}_150.zip")
    return rasterio.open(f"/vsizip/{path:}/{file:}_20101117_gmted_med150.tif")


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
    return 0.11 + 0.67 * numpy.exp(-(slope + 2.0) ** 2 / 1800.)


def all_pairwise_distances(
        elevation: numpy.array,
        transform: rasterio.Affine,
        terrain_coefficients: numpy.array,
):
    d_n, d_e, d_ne = [], [], []
    for y in range(1, len(elevation) + 1):
        (lon0, lat0) = transform * (0, y)
        (lon1, lat1) = transform * (1, y - 1)

        d = GEODESIC.inverse((lon0, lat0), [
            (lon0, lat1), (lon1, lat0), (lon1, lat1)])
        d_n.append(d[0, 0])
        d_e.append(d[1, 0])
        d_ne.append(d[2, 0])
    distance_to_north = numpy.array(d_n)[:-1]
    slope_to_north = 100 * (elevation[1:, :] - elevation[:-1, :]) / distance_to_north[:, None]
    tc_to_north = (terrain_coefficients[1:, :] + terrain_coefficients[:-1, :]) / 2
    north = numpy.full(elevation.shape, numpy.nan)
    north[1:, :] = distance_to_north[:, None]
    north[1:, :] /= (navigation_speed(slope_to_north) * tc_to_north)
    south = numpy.full(elevation.shape, numpy.nan)
    south[:-1, :] = distance_to_north[:, None]
    south[:-1, :] /= (navigation_speed(-slope_to_north) * tc_to_north)
    del distance_to_north, slope_to_north, tc_to_north

    distance_to_east = numpy.array(d_e)
    slope_to_east = 100 * (elevation[:, 1:] - elevation[:, :-1]) / distance_to_east[:, None]
    tc_to_east = (terrain_coefficients[:, 1:] + terrain_coefficients[:, :-1]) / 2
    east = numpy.full(elevation.shape, numpy.nan)
    east[:, :-1] = distance_to_east[:, None]
    east[:, :-1] /= (navigation_speed(slope_to_east) * tc_to_east)
    west = numpy.full(elevation.shape, numpy.nan)
    west[:, 1:] = distance_to_east[:, None]
    west[:, 1:] /= (navigation_speed(-slope_to_east) * tc_to_east)
    del distance_to_east, slope_to_east, tc_to_east

    distance_to_northeast = numpy.array(d_ne)[:-1]
    slope_to_northeast = 100 * (elevation[1:, 1:] - elevation[:-1, :-1]) / distance_to_northeast[:, None]
    tc_to_northeast = (terrain_coefficients[1:, 1:] + terrain_coefficients[:-1, :-1]) / 2
    northeast = numpy.full(elevation.shape, numpy.nan)
    northeast[1:, :-1] = distance_to_northeast[:, None]
    northeast[1:, :-1] /= (navigation_speed(slope_to_northeast) * tc_to_northeast)
    southwest = numpy.full(elevation.shape, numpy.nan)
    southwest[:-1, 1:] = distance_to_northeast[:, None]
    southwest[:-1, 1:] /= (navigation_speed(-slope_to_northeast) * tc_to_northeast)
    del distance_to_northeast, slope_to_northeast, tc_to_northeast
    distance_to_northwest = numpy.array(d_ne)[:-1]
    slope_to_northwest = 100 * (elevation[1:, :-1] - elevation[:-1, 1:]) / distance_to_northwest[:, None]

    tc_to_northwest = (terrain_coefficients[1:, :-1] + terrain_coefficients[:-1, 1:]) / 2

    southeast = numpy.full(elevation.shape, numpy.nan)
    southeast[:-1, :-1] = distance_to_northwest[:, None]
    southeast[:-1, :-1] /= (navigation_speed(-slope_to_northwest) * tc_to_northwest)
    northwest = numpy.full(elevation.shape, numpy.nan)
    northwest[1:, 1:] = distance_to_northwest[:, None]
    northwest[1:, 1:] /= (navigation_speed(slope_to_northwest) * tc_to_northwest)
    del distance_to_northwest, slope_to_northwest, tc_to_northwest

    return {
        (-1, 0): north,
        (-1, 1): northeast,
        (0, 1): east,
        (1, 1): southeast,
        (1, 0): south,
        (1, -1): southwest,
        (0, -1): west,
        (-1, -1): northwest
    }


def travel_time_raster(
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
    bf: int = 500 # Buffer
    elevation_file = gmted_tile_from_geocoordinates(lon, lat)
    profile = elevation_file.profile
    m_trafo = elevation_file.transform
    height, width = elevation_file.shape
    elevation = numpy.full((bf + height + bf, bf + width + bf), -100, int)
    ecoregions = numpy.full((bf + height + bf, bf + width + bf), 999, int)
    for delta_lon, target_cols, source_cols in [
            (-30, slice(0, bf), slice(width-bf, width)),
            (0, slice(bf, width+bf), slice(0, width)),
            (+30, slice(width+bf, width+2*bf), slice(0, bf))]:
        for delta_lat, target_rows, source_rows in [
                (+20, slice(0, bf), slice(height-bf, height)),
                (0, slice(bf, height+bf), slice(0, height)),
                (-20, slice(height+bf, height+2*bf), slice(0, bf))]:
            try:
                elevation[target_rows, target_cols] = (
                    gmted_tile_from_geocoordinates(
                        lon + delta_lon, lat + delta_lat)).read(1)[
                            source_rows, source_cols]
            except rasterio.RasterioIOError:
                pass

    print("Computing hex extents…")
    transform = rasterio.Affine(
        m_trafo.a, 0, m_trafo.c - bf * m_trafo.a,
        0, m_trafo.e, m_trafo.f - bf * m_trafo.e)
    profile["width"] = width + 2 * bf
    profile["height"] = height + 2 * bf
    profile["transform"] = transform
    profile["dtype"] = rasterio.float64
    profile["count"] = 8
    print(profile)

    print("Computing terrain coefficients…")
    terrain_coefficient_raster = TC[ecoregions]

    print("Computing distances on the grid…")
    distance_by_direction = all_pairwise_distances(
        elevation, transform, terrain_coefficient_raster)

    fname = '{:}{:}.tif'.format(lon, lat)
    with rasterio.open(fname, 'w', **profile) as dst:
        for i, band in enumerate(distance_by_direction.values(), 1):
            dst.write(band.astype(rasterio.float64), i)

    # At the end of the ``with rasterio.Env()`` block, context
    # manager exits and all drivers are de-registered.
    return rasterio.open(fname)

class RiverNetwork:
    cache = None
    @classmethod
    def reaches(cls):
        """
        >>> rivers = RiverNetwork.reaches()
        >>> rivers.fields
        [('DeletionFlag', 'C', 1, 0), ['Reach_ID', 'N', 10, 0], ['Next_down', 'N', 10, 0], ['Length_km', 'F', 13, 11], ['Log_Q_avg', 'F', 13, 11], ['Log_Q_var', 'F', 13, 11], ['Class_hydr', 'N', 10, 0], ['Temp_min', 'F', 13, 11], ['CMI_indx', 'F', 13, 11], ['Log_elev', 'F', 13, 11], ['Class_phys', 'N', 10, 0], ['Lake_wet', 'N', 10, 0], ['Stream_pow', 'F', 13, 11], ['Class_geom', 'N', 10, 0], ['Reach_type', 'N', 10, 0], ['Kmeans_30', 'N', 10, 0]]
        >>> rivers.numRecords
        ...

        """
        if cls.cache is not None:
            return cls.cache
        zipshape = zipfile.ZipFile(
            (Path(__file__).parent /
             "../rivers/GloRiC_v10_shapefile.zip").open("rb"))
        shape = shapefile.Reader(
            shp=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.shp"),
            shx=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.shx"),
            dbf=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.dbf"),
            encoding='utf-8'
        )
        cls.cache = shape
        return shape

    def __init__(self, mask: t.Optional[sgeom.Polygon] = None):
        """

        >>> eco = RiverNetwork()
        >>> eco.record(508)[1]
        'Tocantins/Pindare moist forests'

        """
        self.shp = self.reaches()

RIVERS = RiverNetwork()

oceans_zip = zipfile.ZipFile(
            (Path(__file__).parent /
             "../naturalearth/ne_10m_ocean.zip").open("rb"))
OCEANS = shapefile.Reader(
            shp=oceans_zip.open("ne_10m_ocean.shp"),
            shx=oceans_zip.open("ne_10m_ocean.shx"),
            dbf=oceans_zip.open("ne_10m_ocean.dbf"),
            encoding='utf-8'
        )


# Somewhere, I found speeds of 4.5 knots for kayak cruising. That's 8.334 km/h, but the database stores data in seconds.
KAYAK_SPEED = 8.334 / 3600

def estimate_flow_speed(discharge, slope):
    """Estimate the flow speed, in km/s from discharge and slope

    This is a very rough estimate, following [@schulze2005simulating]. They
    suggest to at least estimate the widely varying river roughness n, but we
    cannot do that, so we use their reported mean value of 0.044.

    W = 2.71 · Q^0.557
    D = 0.349 · Q^0.341
    R = D · W / (2 D + W)
    n = 0.044
    v = 1/n · R^2/3 · S^1/2

    """
    n = 0.044
    w = 2.71 * discharge ** 0.557
    d = 0.349 * discharge ** 0.341
    r = d * w / (2 * d + w)
    v = 1/n * r ** (2/3) * slope ** 0.5
    return v / 1000


def navigable_water_mask(dist):
    # Is this reach navigable by Kayak? From
    # [@rood2006instream,@zinke2018comparing] it seems that reaches with a
    # flow lower than 5m³/s are not navigable even by professional extreme
    # sport athletes, and generally even that number seems to be an outlier
    # with opinions starting at 8m³/s, so we take that as the cutoff.
    raster = rasterio.features.rasterize(
        (shp.shape for shp in RIVERS.shp.iterShapeRecords()
         if shp.record[3] > 0.9030899869919434), # log(8.)/log(10.)
        out_shape = dist.shape,
        transform = dist.transform).astype(bool)
    raster |= rasterio.features.rasterize(
        OCEANS.iterShapes(),
        out_shape = dist.shape,
        transform = dist.transform).astype(bool)

    d_n, d_e, d_ne = [], [], []
    for y in range(1, dist.shape[0] + 1):
        (lon0, lat0) = dist.transform * (0, y)
        (lon1, lat1) = dist.transform * (1, y - 1)

        d = GEODESIC.inverse((lon0, lat0), [
            (lon0, lat1), (lon1, lat0), (lon1, lat1)])
        d_n.append(d[0, 0])
        d_e.append(d[1, 0])
        d_ne.append(d[2, 0])

    with rasterio.open("adj.tif", 'w', dist.profile) as dst:
        cell_distance = numpy.ones(dist.shape) * numpy.array(d_n)[:-1, None] / KAYAK_SPEED
        mask = raster[1:, :] & raster[:-1, :]
        with_rivers = dist.read(1)
        dist[1:, :][mask] = cell_distance[mask]
        dst.write(band.astype(rasterio.float64), i)
        south = dist.read(1)
        south[:-1, :][mask] = cell_distance[mask]


    # At the end of the ``with rasterio.Env()`` block, context
    # manager exits and all drivers are de-registered.
    return rasterio.open(fname)



for lon in [-75, -45]:
    for lat in [-20, 0, 20]:
        navigable_water_mask(gmted_tile_from_geocoordinates(lon, lat))
        travel_time = travel_time_raster(lon, lat)
        
