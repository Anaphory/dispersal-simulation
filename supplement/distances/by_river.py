import zipfile
import typing as t
from pathlib import Path

import sqlalchemy
from more_itertools import windowed

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep

from h3 import h3

from database import db

class RiverNetwork:
    cache = None
    @classmethod
    def reaches(cls):
        """
        >>> rivers = RiverNetwork.reaches()
        >>> rivers.numRecords
        847
        >>> rivers.fields
        [('DeletionFlag', 'C', 1, 0), ['Reach_ID', 'N', 10, 0], ['Next_down', 'N', 10, 0], ['Length_km', 'F', 13, 11], ['Log_Q_avg', 'F', 13, 11], ['Log_Q_var', 'F', 13, 11], ['Class_hydr', 'N', 10, 0], ['Temp_min', 'F', 13, 11], ['CMI_indx', 'F', 13, 11], ['Log_elev', 'F', 13, 11], ['Class_phys', 'N', 10, 0], ['Lake_wet', 'N', 10, 0], ['Stream_pow', 'F', 13, 11], ['Class_geom', 'N', 10, 0], ['Reach_type', 'N', 10, 0], ['Kmeans_30', 'N', 10, 0]]

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



engine, tables = db('sqlite:///rivers.sqlite')
t_hex = tables["hex"]
t_reach = tables["reach"]
t_flows = tables["flows"]

downstream_from_navigable = set()

for r, reach in enumerate(RIVERS.shp.iterShapeRecords()):
    points: t.List[t.Tuple[float, float]] = reach.shape.points
    data = reach.record
    reach_id = int(data[0])

    # Is this reach navigable by Kayak? From
    # [@rood2006instream,@zinke2018comparing] it seems that reaches with a flow
    # lower than 5m³/s are not navigable even by professional extreme sport
    # athletes, and generally even that number seems to be an outlier with
    # options starting at 8m³/s, so we take that as the cutoff.
    #
    # [@zinke2018comparing] further plots wild water kayaking run slopes vs.
    # difficulty. All of these are below 10%, so we assume that reaches above
    # 10% are not navigable. Gradient is not directly available in the data,
    # but the stream power is directly proportional to the product of discharge
    # and gradient, so we can reverse-engineer it:
    # Stream Power [kg m/s³]
    # = Water Density [kg/m³] * gravity [m/s²] * discharge [m³/s] * slope [m/m]
    # so
    # slope = stream power / discharge / (1000 * 9.81) > 10% = 0.1
    if reach_id not in downstream_from_navigable or (
            data[3] < 0.9030899869919434 or data[11] / 10 ** data[3] > 981.0):
        continue
    if data[1] < reach_id:
        print(data[1])
    downstream_from_navigable.add(data[1])
    downstream_from_navigable.discard(reach_id)
    print(data)
    try:
        engine.execute(t_reach.insert({"id": reach_id, "index": r}))
    except sqlalchemy.exc.IntegrityError:
        engine.execute(
            t_reach.update().where(
                t_reach.c.id == reach_id).values(
                    {"index": r}))
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
    for hex in hexes:
        try:
            engine.execute(t_hex.insert({"hexbin": hex}))
        except sqlalchemy.exc.IntegrityError:
            pass
        engine.execute(t_flows.insert({"reach": reach_id, "hexbin": hex}))
