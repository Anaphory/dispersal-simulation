import zipfile
import typing as t
from pathlib import Path

import sqlalchemy
from more_itertools import windowed

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep
from shapely.ops import unary_union
import cartopy.geodesic as geodesic
import cartopy.io.shapereader as shpreader

from h3 import h3

from database import db

from by_river import KAYAK_SPEED # in km/s

land_shp_fname = shpreader.natural_earth(
    resolution='50m', category='physical', name='land')

land_geom = unary_union([
    record.geometry
    for record in shpreader.Reader(land_shp_fname).records()
    if record.attributes.get('featurecla') != "Null island"])

LAND = prep(land_geom)

DEFINITELY_INLAND = prep(land_geom.buffer(-1/60., resolution = 4)) # 1 arc minute – 2 km near the equator.
BUFFER_NOT_SEA = prep(land_geom.buffer(1/60., resolution = 4))

GEODESIC = geodesic.Geodesic()

def is_landlocked(xy):
    return DEFINITELY_INLAND.contains(sgeom.Point(*xy))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("--start")
    args = parser.parse_args()
    engine, tables = db(args.database)
    t_hex = tables["hex"]
    t_dist = tables["dist"]

    failed = set()
    query = sqlalchemy.select([t_hex.c.hexbin, t_hex.c.vlongitude, t_hex.c.vlatitude])
    if args.start:
        query = query.where(t_hex.c.hexbin >= args.start)
    for hexbin0, lon0, lat0 in engine.execute(query).fetchall():
        try:
            if is_landlocked((lon0, lat0)):
                continue
        except TypeError:
            failed.add(hexbin0)
            continue
        print("Working on {:}…".format(hexbin0))
        for hexbin1, lon1, lat1 in engine.execute(
            sqlalchemy.select([t_hex.c.hexbin, t_hex.c.vlongitude, t_hex.c.vlatitude])
                .where(t_hex.c.vlongitude > lon0 - 2.5)
                .where(t_hex.c.vlongitude < lon0 + 2.5)
                .where(t_hex.c.vlatitude > lat0 - 2.5)
                .where(t_hex.c.vlatitude < lat0 + 2.5)
        ).fetchall():
            d = GEODESIC.inverse((lon0, lat0), (lon1, lat1))[0, 0]
            if is_landlocked((lon1, lat1)):
                continue
            s = 8
            if d > 300_000: # meters
                continue
            if DEFINITELY_INLAND.intersects(sgeom.LineString([(lon0, lat0), (lon1, lat1)])):
                continue
            if BUFFER_NOT_SEA.intersects(sgeom.LineString([(lon0, lat0), (lon1, lat1)])):
                print("Path from {:} to {:} might cross land.".format(hexbin0, hexbin1))
                s = 7
            t = d / 1000. / KAYAK_SPEED
            if t > 8. * 60. * 24.:
                t *= 3
            try:
                engine.execute(t_dist.insert(
                    {"hexbin1": hexbin0, "hexbin2": hexbin1, "source": s, "distance": t, "flat_distance": d}))
            except sqlalchemy.exc.IntegrityError:
                pass
        print(hexbin0, "done")
    print(failed)
