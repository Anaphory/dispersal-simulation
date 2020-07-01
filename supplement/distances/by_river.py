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

RESOLUTION = 5


def neighbors(hex: h3.H3Index) -> t.Set[h3.H3Index]:
    this, neighbors = h3.k_ring_distances(hex, 1)
    return neighbors


def intersection(xa, ya, xb, yb, xl, yl, xr, yr):
    """

    >>> intersection(0, 1, 3, 1, 2, 3, 2, 0)
    (0.6666666666666666, 0.3333333333333333)
    """

    if xl != xr:
        t = (xl * (ya - yb) + xa * (yb - yl) + xb * (yl - ya)) / (
            (xl - xr) * (ya - yb) + xa * (yr - yl) + xb * (yl - yr))
    else:
        t = (yl * (xa - xb) + ya * (xb - xl) + yb * (xl - xa)) / (
            (yl - yr) * (xa - xb) + ya * (xr - xl) + yb * (xl - xr))

    if xb != xa:
        u = (xb * (yl - yr) + xl * (yr - yb) + xr * (yb - yl)) / (
            (xb - xa) * (yl - yr) + xl * (ya - yb) + xr * (yb - ya))
    else:
        u = (yb * (xl - xr) + yl * (xr - xb) + yr * (xb - xl)) / (
            (yb - ya) * (xl - xr) + yl * (xa - xb) + yr * (xb - xa))

    return t, u


def find_hexes_crossed(
        x0: float, y0: float, x1: float, y1: float,
        length_start: float, length_end:float
) -> t.List[t.Tuple[h3.H3Index, float, float]]:
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
        return [(start_hex, length_start, length_end)]
    elif end_hex in neighbors(start_hex):
        # Check whether the line segment runs through the left neighbor,
        # through the right neighbor, or through the shared boundary of start
        # and end.
        nr, nl = neighbors(end_hex) & neighbors(start_hex)
        ya, xa = h3.h3_to_geo(start_hex)
        yb, xb = h3.h3_to_geo(end_hex)
        yr, xr = h3.h3_to_geo(nr)
        yl, xl = h3.h3_to_geo(nl)

        # strange things can happen when wrapping the date line.
        t, u = intersection(xa, ya, xb, yb, xl, yl, xr, yr)
        assert 1./6. < t < 5./6.

        if u < 1e-10:
            # The line practially runs completely on the side of the end hex
            return [(end_hex, length_start, length_end)]
        if u > 1-1e-10:
            # The line practially runs completely on the side of the start hex
            return [(start_hex, length_start, length_end)]

        crossing_point = length_start + u * (length_end - length_start)

        if t < 1./3.:
            x_m, y_m = xa + u * (xb - xa), yb + u * (yb - ya)
            return find_hexes_crossed(x0, y0, x_m, y_m, length_start, crossing_point) + find_hexes_crossed(x_m, y_m, x1, y1, crossing_point, length_end)
        elif t > 2./3.:
            x_m, y_m = xa + u * (xb - xa), yb + u * (yb - ya)
            return find_hexes_crossed(x0, y0, x_m, y_m, length_start, crossing_point) + find_hexes_crossed(x_m, y_m, x1, y1, crossing_point, length_end)
        else:
            return [(start_hex, length_start, crossing_point), (end_hex, crossing_point, length_end)]
    else:
        xmid = (x0 + x1) / 2.
        ymid = (y0 + y1) / 2.
        length_mid = length_start + 0.5 * (length_end - length_start)
        return find_hexes_crossed(x0, y0, xmid, ymid, length_start, length_mid) + find_hexes_crossed(xmid, ymid, x1, y1, length_mid, length_end)



engine, tables = db('sqlite:///run.sqlite')
t_hex = tables["hex"]
t_reach = tables["reach"]
t_flows = tables["flows"]

downstream_from_navigable = set()


def estimate_flow_speed(discharge, slope):
    """Estimate the flow speed, in km/h from discharge and slope

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
    return v * 3.6

if __name__ == "__main__":
    for r, reach in enumerate(RIVERS.shp.iterShapeRecords()):
        data = reach.record
        reach_id = int(data[0])
        if reach_id < 60000000:
            # The Americas have SA 6, NA 7, American Arctic 8, Greenland 9
            continue

        points: t.List[t.Tuple[float, float]] = reach.shape.points

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
        if reach_id not in downstream_from_navigable and (
                data[3] < 0.9030899869919434 or data[11] / 10 ** data[3] > 981.0):
            continue
        if data[1] < reach_id:
            print(data[1])
        downstream_from_navigable.add(data[1])
        downstream_from_navigable.discard(reach_id)
        v = estimate_flow_speed(discharge = 10 ** data[3], slope = data[11] / 10 ** data[3] / (9810))
        try:
            engine.execute(t_reach.insert(
                {"reach_id": reach_id, "record_index": r,
                "next_down": data[1],
                # Somewhere, I found speeds of 4.5 knots for kayak cruising.
                "travel_downstream": data[2] / (8.334 + v),
                "travel_upstream": data[2] / max(8.334 - v, 0.3)}))
        except sqlalchemy.exc.IntegrityError:
            engine.execute(t_reach.update(
                { "record_index": r,
                  "next_down": data[1],
                  # Somewhere, I found speeds of 4.5 knots for kayak cruising.
                  "travel_downstream": data[2] / (8.334 + v),
                  "travel_upstream": data[2] / max(8.334 - v, 0.3)}).where(
                      "reach_id" == reach_id,
                  ))
            continue

        hexes: t.Dict[h3.H3Index, t.Tuple[float, float]] = {}
        length = 0.0
        for start, end in windowed(points, 2):
            if start is None:
                continue
            else:
                x0, y0 = start
            if end is None:
                continue
            else:
                x1, y1 = end
            length_after = length + ((x0 - x1) ** 2 + (y0 - y1) ** 2) ** 0.5
            try:
                for h, st, ed in find_hexes_crossed(x0, y0, x1, y1, length, length_after):
                    if h in hexes:
                        if hexes[h][0] > st:
                            hexes[h] = (st, hexes[h][1])
                        if hexes[h][1] < ed:
                            hexes[h] = (hexes[h][0], ed)
                    else:
                        hexes[h] = st, ed
            except AssertionError:
                # Probably near poles or date line.
                continue
            length = length_after
    
        for hex, (st, ed) in hexes.items():
            try:
                engine.execute(t_hex.insert({"hexbin": hex}))
            except sqlalchemy.exc.IntegrityError:
                pass
            engine.execute(t_flows.insert(
                {"reach_id": reach_id, "hexbin": hex, "start": st/length, "end": ed/length}))
    print(reach_id)
