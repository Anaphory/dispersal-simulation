from tqdm import tqdm
import zipfile
import typing as t
from pathlib import Path
import numpy

import sqlalchemy
from more_itertools import windowed

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep
from shapely.ops import unary_union
import cartopy.geodesic as geodesic
import cartopy.io.shapereader as shpreader

from h3.api import basic_int as h3

from database import db

from by_river import KAYAK_SPEED  # in m/s

GEODESIC: geodesic.Geodesic = geodesic.Geodesic()

DATABASE, TABLES = db()

def distance_by_sea(definitely_inland):
    query = sqlalchemy.select(
        [TABLES["nodes"].c.node_id, TABLES["nodes"].c.longitude, TABLES["nodes"].c.latitude]
    ).where(TABLES["nodes"].c.coastal)
    for node0, lon0, lat0 in tqdm(DATABASE.execute(query).fetchall()):
        with DATABASE.begin() as conn:
            for node1, lon1, lat1 in DATABASE.execute(
                query.where(lat0 - TABLES["nodes"].c.latitude > - 3.0)
                .where(lat0 - TABLES["nodes"].c.latitude < 3.0)
                .where(
                    (lon0 - TABLES["nodes"].c.longitude) * numpy.cos(lat0 * numpy.pi / 180) < 3.0
                )
                .where(
                    (lon0 - TABLES["nodes"].c.longitude) * numpy.cos(lat0 * numpy.pi / 180) > -3.0
                )
            ).fetchall():
                d = GEODESIC.inverse((lon0, lat0), (lon1, lat1))[0, 0]
                if d > 300_000:  # meters
                    continue
                if definitely_inland.intersects(
                    sgeom.LineString([(lon0, lat0), (lon1, lat1)])
                ):
                    continue
                t = d / KAYAK_SPEED
                conn.execute(
                    TABLES["edges"].insert(
                        {
                            "node1": node0,
                            "node2": node1,
                            "source": "sea",
                            "travel_time": t,
                            "flat_distance": d,
                        }
                    )
                )
