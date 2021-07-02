import typing as t
from pathlib import Path

import rasterio

Tile = t.Tuple[t.Literal["N", "S"], int, t.Literal["E", "W"], int]
RowCol = t.Tuple[int, int]


def gmted_tile(tile: Tile) -> rasterio.DatasetReader:
    """Path to the GMTED tile for given geocoordinates.

    The return type is explicitly a `str`, not a `pathlib.Path`, because it
    contains a GDAL virtual file system component.

    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010N30E000/
    >>> d1 = gmted_tile(tile_from_geocoordinates(15, 40))
    >>> d1.name[-78:]
    'elevation/GMTED2010/GMTED2010N30E000_150.zip/30n000e_20101117_gmted_med150.tif'

    https://earthexplorer.usgs.gov/metadata/4584/GMTED2010S30W090/
    >>> d2 = gmted_tile(tile_from_geocoordinates(-75, -20))
    >>> d2.name[-78:]
    'elevation/GMTED2010/GMTED2010S30W090_150.zip/30s090w_20101117_gmted_med150.tif'

    """
    ns, lat, ew, lon = tile
    file = "{:02d}{:1s}{:03d}{:1s}".format(lat, ns.lower(), lon, ew.lower())
    tile = "GMTED2010{:1s}{:02d}{:1s}{:03d}".format(ns.upper(), lat, ew.upper(), lon)
    path = (
        Path(__file__).absolute().parent.parent / f"elevation/GMTED2010/{tile:}_150.zip"
    )
    return rasterio.open(f"/vsizip/{path:}/{file:}_20101117_gmted_med150.tif")


def ecoregion_tile(tile: Tile) -> rasterio.DatasetReader:
    ns, lat0, ew, lon0 = tile
    ecoregions_path_t = (
        "../ecoregions/ECOREGIONS-{0:02d}{1}{2:03d}{3}_20101117_gmted_med150.tif"
    )
    return rasterio.open(ecoregions_path_t.format(lat0, ns.lower(), lon0, ew.lower()))


def tile_from_geocoordinates(lon: float, lat: float) -> Tile:
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


def boundingbox_from_tile(tile: Tile):
    west = (-1 if tile[2] == "W" else 1) * tile[3]
    south = (-1 if tile[0] == "S" else 1) * tile[1]
    return (west, south, west + 30, south + 20)
