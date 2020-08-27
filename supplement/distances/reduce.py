"""
reduce – Remove uninhabitable nodes while keeping the graph structure
"""

import json
import numpy
import zipfile
import collections
import typing as t
from pathlib import Path
from itertools import count
from dataclasses import dataclass
from heapq import heappush as push, heappop as pop

import numpy
import sqlalchemy
from sqlalchemy.sql import func

import matplotlib.pyplot as plt
import matplotlib.collections
import matplotlib.transforms as mtransforms

import cartopy.geodesic as geodesic
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

import shapefile
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

import rasterio

from h3 import h3
from h3.h3 import H3Index as Index

from database import db
from ecoregions import ECOREGIONS, TC, RESOLUTION
from distances import lonlat

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("--tile", action="append", default=[], type=lonlat)
    args = parser.parse_args()
    engine, tables = db(args.database)
    t_hex = tables["hex"]
    t_dist = tables["dist"]
    t_dist1 = t_dist.alias("dist1")
    t_dist2 = t_dist.alias("dist2")
    t_eco = tables["eco"]

    full = sqlalchemy.select(
        [sqlalchemy.text("eco.hexbin AS habitat")]
    ).where(
        t_eco.c.ecoregion != 999
    ).group_by(
        t_eco.c.hexbin
    )
    empty = sqlalchemy.select([t_hex.c.hexbin]).select_from(
        t_hex.join(full,
                   onclause=sqlalchemy.text("habitat is NULL"),
                   isouter=True))
    for i, in engine.execute(empty).fetchall():
        print(i)
        query = sqlalchemy.select([
            t_dist.c.hexbin1, t_dist1.c.hexbin2,
            t_dist.c.distance + t_dist1.c.distance,
            t_dist.c.source * t_dist1.c.source]
        ).select_from(
            t_dist.join(
                t_dist1,
                onclause=(
                    (t_dist.c.hexbin2 == i) &
                    (t_dist1.c.hexbin1 == i))
            ).join(
                t_dist2,
                onclause=((t_dist2.c.hexbin1 == t_dist.c.hexbin1) &
                          (t_dist2.c.hexbin2 == t_dist1.c.hexbin2)),
                isouter=True
            )
        ).where(
            ((t_dist2.c.hexbin1 == None) &
             (t_dist2.c.hexbin2 == None)) |
            (t_dist2.c.distance > t_dist.c.distance + t_dist1.c.distance) &
            (t_dist.c.distance + t_dist1.c.distance < 3*8*60*60)
        )
        shortcuts = t_dist.insert().from_select(
            [t_dist.c.hexbin1, t_dist.c.hexbin2,
             t_dist.c.distance, t_dist.c.source],
            query).prefix_with("OR REPLACE")
        with engine.begin() as conn:
            conn.execute(shortcuts)
            conn.execute(t_eco.delete(whereclause=t_eco.c.hexbin == i))
            conn.execute(t_dist.delete(whereclause=t_dist.c.hexbin1 == i))
            conn.execute(t_dist.delete(whereclause=t_dist.c.hexbin2 == i))
            conn.execute(t_hex.delete(whereclause=t_hex.c.hexbin == i))
