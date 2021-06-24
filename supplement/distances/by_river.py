import zipfile
import typing as t
from pathlib import Path

from tqdm import tqdm

import sqlalchemy
from sqlalchemy.dialects.sqlite import insert

from more_itertools import windowed

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep

from h3.api import basic_int as h3

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
            (Path(__file__).parent / "../rivers/GloRiC_v10_shapefile.zip").open("rb")
        )
        shape = shapefile.Reader(
            shp=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.shp"),
            shx=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.shx"),
            dbf=zipshape.open("GloRiC_v10_shapefile/GloRiC_v10.dbf"),
            encoding="utf-8",
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


RESOLUTION = 5


def neighbors(hex) -> t.Set:
    this, neighbors = h3.k_ring_distances(hex, 1)
    return neighbors


def intersection(xa, ya, xb, yb, xl, yl, xr, yr):
    """

    >>> intersection(0, 1, 3, 1, 2, 3, 2, 0)
    (0.6666666666666666, 0.3333333333333333)
    """

    if xl != xr:
        t = (xl * (ya - yb) + xa * (yb - yl) + xb * (yl - ya)) / (
            (xl - xr) * (ya - yb) + xa * (yr - yl) + xb * (yl - yr)
        )
    else:
        t = (yl * (xa - xb) + ya * (xb - xl) + yb * (xl - xa)) / (
            (yl - yr) * (xa - xb) + ya * (xr - xl) + yb * (xl - xr)
        )

    if xb != xa:
        u = (xb * (yl - yr) + xl * (yr - yb) + xr * (yb - yl)) / (
            (xb - xa) * (yl - yr) + xl * (ya - yb) + xr * (yb - ya)
        )
    else:
        u = (yb * (xl - xr) + yl * (xr - xb) + yr * (xb - xl)) / (
            (yb - ya) * (xl - xr) + yl * (xa - xb) + yr * (xb - xa)
        )

    return t, u


def find_hexes_crossed(
    x0: float, y0: float, x1: float, y1: float, length_start: float, length_end: float
) -> t.List[t.Tuple]:
    """Find all hexes that a line segment crosses.

    By bisecting the line segment from (x0, y0) to (x1, y1), find all h3 hexes
    that the line crosses.

    >>> find_hexes_crossed(-10, -10, -9, -9)
    ...

    """
    start_hex = h3.geo_to_h3(y0, x0, RESOLUTION)
    end_hex = h3.geo_to_h3(y1, x1, RESOLUTION)
    if start_hex == end_hex:
        # Hexes are convex, so the line segment must run entirely within.
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
        assert 1.0 / 6.0 < t < 5.0 / 6.0

        if u < 1e-10:
            # The line practially runs completely on the side of the end hex
            return [(end_hex, length_start, length_end)]
        if u > 1 - 1e-10:
            # The line practially runs completely on the side of the start hex
            return [(start_hex, length_start, length_end)]

        crossing_point = length_start + u * (length_end - length_start)

        if t < 1.0 / 3.0:
            x_m, y_m = xa + u * (xb - xa), yb + u * (yb - ya)
            return find_hexes_crossed(
                x0, y0, x_m, y_m, length_start, crossing_point
            ) + find_hexes_crossed(x_m, y_m, x1, y1, crossing_point, length_end)
        elif t > 2.0 / 3.0:
            x_m, y_m = xa + u * (xb - xa), yb + u * (yb - ya)
            return find_hexes_crossed(
                x0, y0, x_m, y_m, length_start, crossing_point
            ) + find_hexes_crossed(x_m, y_m, x1, y1, crossing_point, length_end)
        else:
            return [
                (start_hex, length_start, crossing_point),
                (end_hex, crossing_point, length_end),
            ]
    else:
        xmid = (x0 + x1) / 2.0
        ymid = (y0 + y1) / 2.0
        length_mid = length_start + 0.5 * (length_end - length_start)
        return find_hexes_crossed(
            x0, y0, xmid, ymid, length_start, length_mid
        ) + find_hexes_crossed(xmid, ymid, x1, y1, length_mid, length_end)


# 3 knots is about 1.5433333 m/s
KAYAK_SPEED = 1.5433333


def estimate_flow_speed(discharge, slope):
    """Estimate the flow speed, in m/s from discharge (m³/s) and slope (m/m)

    This is a very rough estimate, following [@schulze2005simulating]. They
    suggest to at least estimate the widely varying river roughness n, but we
    cannot do that, so we use their reported mean value of 0.044.

    W = 2.71 · Q^0.557
    D = 0.349 · Q^0.341
    R = D · W / (2 D + W)
    n = 0.044

    # Manning-Strickler formula
    v = 1/n · R^2/3 · S^1/2

    """
    n = 0.044
    w = 2.71 * discharge ** 0.557
    d = 0.349 * discharge ** 0.341
    r = d * w / (2 * d + w)
    # Manning-Strickler formula
    v = 1 / n * r ** (2 / 3) * slope ** 0.5
    return v


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("--start", type=int)
    args = parser.parse_args()
    engine, tables = db(args.database)
    t_node = tables["nodes"]
    t_dist = tables["edges"]

    RIVERS = RiverNetwork()

    for r, reach in tqdm(enumerate(RIVERS.shp.iterShapeRecords())):
        data = reach.record
        reach_id = int(data[0])
        if reach_id < 60000000:
            # The Americas have SA 6, NA 7, American Arctic 8, Greenland 9
            continue

        points: t.List[t.Tuple[float, float]] = reach.shape.points

        # Is this reach navigable by Kayak? From
        # [@rood2006instream,@zinke2018comparing] it seems that reaches with a
        # flow lower than 5m³/s are not navigable even by professional extreme
        # sport athletes, and generally even that number seems to be an outlier
        # with opinions starting at 8m³/s, so we take that as the cutoff.
        #
        # [@zinke2018comparing] further plots wild water kayaking run slopes
        # vs. difficulty. All of these are below 10%, so we assume that reaches
        # above 10% are not navigable. Gradient is not directly available in
        # the data, but the stream power is directly proportional to the
        # product of discharge and gradient, so we can reverse-engineer it:
        # Stream Power [kg m/s³]
        # = Water Density [kg/m³] * gravity [m/s²] * discharge [m³/s] * slope [m/m]
        # so
        if (
            data[3] < 0.9030899869919434
            or data[11] / (10 ** data[3]) > 981.0  # log(8.)/log(10.)
        ):
            continue

        with engine.begin() as conn:
            # slope = stream power / discharge / (1000 * 9.81) > 10% = 0.1
            v = estimate_flow_speed(
                discharge=10 ** data[3], slope=data[11] / 10 ** data[3] / (9810)
            )
            print(data[0], KAYAK_SPEED + v, KAYAK_SPEED + v)
            if v > 40.0:
                raise RuntimeError(
                    "Found a reach flowing faster than 40 m/s, something is clearly wrong."
                )

            downstream = data[1]
            conn.execute(
                insert(t_node)
                .values(
                    {
                        "node_id": reach_id,
                        "longitude": points[0][0],
                        "latitude": points[0][1],
                        "coastal": False,
                    },
                )
                .on_conflict_do_nothing()
            )
            if downstream == 0:
                downstream = -reach_id
                coastal = True
            else:
                coastal = False
            conn.execute(
                insert(t_node)
                .values(
                    {
                        "node_id": downstream,
                        "longitude": points[-1][0],
                        "latitude": points[-1][1],
                        "coastal": coastal,
                    },
                )
                .on_conflict_do_nothing()
            )

            conn.execute(
                insert(t_dist)
                .values(
                    {
                        "node1": reach_id,
                        "node2": downstream,
                        "source": "river",
                        "flat_distance": data[2] * 1000,
                        "travel_time": data[2] * 1000 / (KAYAK_SPEED + v),
                    },
                )
                .on_conflict_do_nothing()
            )
            if KAYAK_SPEED > v:
                conn.execute(
                    insert(t_dist)
                    .values(
                        {
                            "node1": downstream,
                            "node2": reach_id,
                            "source": "river",
                            "flat_distance": data[2] * 1000,
                            "travel_time": data[2] * 1000 / (KAYAK_SPEED - v),
                        },
                    )
                    .on_conflict_do_nothing()
                )
