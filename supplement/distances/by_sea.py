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


def compute():
    failed = set()
    query = sqlalchemy.select(
        [t_hex.c.node_id, t_hex.c.longitude, t_hex.c.latitude]
    ).where(t_hex.c.coastal)
    for node0, lon0, lat0 in tqdm(engine.execute(query).fetchall()):
        with engine.begin() as conn:
            for node1, lon1, lat1 in engine.execute(
                sqlalchemy.select(
                    [t_hex.c.node_id, t_hex.c.longitude, t_hex.c.latitude]
                )
                .where(t_hex.c.latitude > lat0 - 2.0)
                .where(t_hex.c.latitude < lat0 + 2.0)
                .where(
                    t_hex.c.longitude > lon0 - 2.0 / numpy.cos(lat0 * numpy.pi / 180)
                )
                .where(
                    t_hex.c.longitude < lon0 + 2.0 / numpy.cos(lat0 * numpy.pi / 180)
                )
                .where(t_hex.c.coastal)
            ).fetchall():
                d = GEODESIC.inverse((lon0, lat0), (lon1, lat1))[0, 0]
                if d > 300_000:  # meters
                    continue
                if DEFINITELY_INLAND.intersects(
                    sgeom.LineString([(lon0, lat0), (lon1, lat1)])
                ):
                    continue
                t = d / KAYAK_SPEED
                conn.execute(
                    t_dist.insert(
                        {
                            "node1": hexbin0,
                            "node2": hexbin1,
                            "source": "sea",
                            "travel_time": t,
                            "flat_distance": d,
                        }
                    )
                )
