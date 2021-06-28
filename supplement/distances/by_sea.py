#!/home/cluster/gkaipi/.pyenv/shims/python
from tqdm import tqdm
import zipfile
import typing as t
from pathlib import Path
import numpy

import sqlalchemy
from sqlalchemy.dialects.sqlite import insert
from more_itertools import windowed

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep
from shapely.ops import unary_union
import cartopy.geodesic as geodesic
import cartopy.io.shapereader as shpreader
from earth import LAND

from h3.api import basic_int as h3

from database import db

from by_river import KAYAK_SPEED  # in m/s

GEODESIC: geodesic.Geodesic = geodesic.Geodesic()

DATABASE, TABLES = db(
    file = "sqlite:///migration-network-add-sea.sqlite",
)

def distance_by_sea(definitely_inland, skip=True):
    query = sqlalchemy.select(
        [
            TABLES["nodes"].c.node_id,
            TABLES["nodes"].c.longitude,
            TABLES["nodes"].c.latitude,
        ]
    ).where(TABLES["nodes"].c.coastal)
    for node0, lon0, lat0 in tqdm(DATABASE.execute(query).fetchall()):
        if (
            skip
            and DATABASE.execute(
                sqlalchemy.select([TABLES["edges"].c.node1]).where(
                    TABLES["edges"].c.node1 == node0,
                    TABLES["edges"].c.travel_time != None,
                    TABLES["edges"].c.source == "sea",
                )
            ).fetchall()
        ):
            continue
        DATABASE.execute(
            insert(TABLES["edges"])
            .values(
                node1=node0,
                node2=node0,
                source="sea",
                travel_time=0.0,
                flat_distance=0.0,
            )
            .on_conflict_do_nothing()
        )
        for node1, lon1, lat1 in DATABASE.execute(
                query.where((lat0 - TABLES["nodes"].c.latitude ) ** 2
                   + (lon0 - TABLES["nodes"].c.longitude) ** 2
                    * numpy.cos(lat0 * numpy.pi / 180)**2
                    < 9)
            ).fetchall():
                d = GEODESIC.inverse((lon0, lat0), (lon1, lat1))[0, 0]
                if d > 300_000:  # meters
                    continue
                if definitely_inland.intersects(
                    sgeom.LineString([(lon0, lat0), (lon1, lat1)])
                ):
                    continue
                t = d / KAYAK_SPEED
                DATABASE.execute(
                    insert(TABLES["edges"])
                    .values(
                        node1=node0,
                        node2=node1,
                        source="sea",
                        travel_time=t,
                        flat_distance=d,
                    )
                    .on_conflict_do_nothing()
                )

distance_by_sea(LAND.buffer(-0.04))
