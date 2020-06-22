import zipfile
import typing as t
from pathlib import Path

from more_itertools import windowed

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep

from h3 import h3


class RiverNetwork:
    cache = None
    @classmethod
    def reaches(cls):
        """
        >>> rivers = RiverNetwork.reaches()
        >>> eco.numRecords
        847
        >>> eco.fields
        [('DeletionFlag', 'C', 1, 0), ['OBJECTID', 'N', 32, 10], ['ECO_NAME', 'C', 150, 0], ['BIOME_NUM', 'N', 32, 10], ['BIOME_NAME', 'C', 254, 0], ['REALM', 'C', 254, 0], ['ECO_BIOME_', 'C', 254, 0], ['NNH', 'N', 11, 0], ['ECO_ID', 'N', 11, 0], ['SHAPE_LENG', 'N', 32, 10], ['SHAPE_AREA', 'N', 32, 10], ['NNH_NAME', 'C', 64, 0], ['COLOR', 'C', 7, 0], ['COLOR_BIO', 'C', 7, 0], ['COLOR_NNH', 'C', 7, 0], ['LICENSE', 'C', 64, 0]]

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

RESOLUTION = 5


def neighbors(hex: h3.H3Index) -> t.Set[h3.H3Index]:
    this, neighbors = h3.k_ring_distances(hex, 1)
    return neighbors


def find_hexes_crossed(
        x0: float, y0: float, x1: float, y1: float
) -> t.Set[h3.H3Index]:
    """Find all hexes that a line segment crosses.

    By bisecting the line segment from (x0, y0) to (x1, y1), find all h3 hexes
    that the line crosses.

    >>> find_hexes_crossed(-10, -10, -9, -9)
    ...

    """
    start_hex: h3.H3Index = h3.geo_to_h3(y0, x0, RESOLUTION)
    end_hex: h3.H3Index = h3.geo_to_h3(y1, x1, RESOLUTION)
    if start_hex == end_hex:
        # Hexes are covex, so the line segment must run entirely within.
        return {start_hex}
    elif end_hex in neighbors(start_hex):
        # Check whether the line segment runs through the left neighbor,
        # through the right neighbor, or through the shared boundary of start
        # and end.
        nr, nl = neighbors(end_hex) & neighbors(start_hex)
        ya, xa = h3.h3_to_geo(start_hex)
        yb, xb = h3.h3_to_geo(end_hex)
        yr, xr = h3.h3_to_geo(nr)
        yl, xl = h3.h3_to_geo(nl)

        if xl != xr:
            t = (xl * (ya - yb) + xa * (yb - yl) + xb * (yl - ya)) / (
                (xl - xr) * (ya - yb) + xa * (yr - yl) + xb * (yl - yr))
        else:
            t = (yl * (xa - xb) + ya * (xb - xl) + yb * (xl - xa)) / (
                (yl - yr) * (xa - xb) + ya * (xr - xl) + yb * (xl - xr))
        if t < 1./3.:
            return {start_hex, end_hex, nl}
        elif t > 2./3.:
            return {start_hex, end_hex, nr}
        else:
            return {start_hex, end_hex}
    else:
        breakpoint()
        xmid = (x0 + x1) / 2.
        ymid = (y0 + y1) / 2.
        return find_hexes_crossed(
            x0, y0, xmid, ymid) | find_hexes_crossed(
                xmid, ymid, x1, y1)


for reach in RIVERS.shp.iterShapeRecords():
    points: t.List[t.Tuple[float, float]] = reach.shape.points
    data = reach.record

    hexes: t.Set[h3.H3Index] = set()
    for start, end in windowed(points, 2):
        if start is None:
            continue
        else:
            x0, y0 = start
        if end is None:
            continue
        else:
            x1, y1 = end
        hexes.update(find_hexes_crossed(x0, y0, x1, y1))

