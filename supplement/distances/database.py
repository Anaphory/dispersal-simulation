import typing as t

import sqlalchemy
from sqlalchemy import event
from sqlalchemy.engine import Engine

from h3.api import basic_int as h3


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(
    dbapi_connection: sqlalchemy.engine.Connectable, connection_record: t.Any
) -> None:
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()


def db(
    file: str = "sqlite:///migration-network.sqlite",
) -> t.Tuple[sqlalchemy.engine.Connectable, t.Dict[str, sqlalchemy.Table]]:
    engine = sqlalchemy.create_engine(file)
    metadata = sqlalchemy.MetaData()

    nodes = sqlalchemy.Table(
        "nodes",
        metadata,
        sqlalchemy.Column("node_id", sqlalchemy.Integer, primary_key=True),
        sqlalchemy.Column("longitude", sqlalchemy.Float),
        sqlalchemy.Column("latitude", sqlalchemy.Float),
        sqlalchemy.Column("h3longitude", sqlalchemy.Float),
        sqlalchemy.Column("h3latitude", sqlalchemy.Float),
        sqlalchemy.Column("coastal", sqlalchemy.Boolean),
    )

    edges = sqlalchemy.Table(
        "edges",
        metadata,
        sqlalchemy.Column(
            "node1", sqlalchemy.Integer, sqlalchemy.ForeignKey(nodes.c.node_id), primary_key=True
        ),
        sqlalchemy.Column(
            "node2", sqlalchemy.Integer, sqlalchemy.ForeignKey(nodes.c.node_id), primary_key=True
        ),
        sqlalchemy.Column("travel_time", sqlalchemy.Float), # in seconds
        sqlalchemy.Column("flat_distance", sqlalchemy.Float), # in meters
        sqlalchemy.Column("source", sqlalchemy.String, primary_key=True),
    )

    ecology = sqlalchemy.Table(
        "ecology",
        metadata,
        sqlalchemy.Column(
            "node", sqlalchemy.Integer, sqlalchemy.ForeignKey(nodes.c.node_id), primary_key=True
        ),
        sqlalchemy.Column("ecoregion", sqlalchemy.Integer(), primary_key=True),
        sqlalchemy.Column("population_capacity", sqlalchemy.Float),
    )

    metadata.create_all(engine)
    return engine, {
        "nodes": nodes,
        "edges": edges,
        "ecology": ecology
    }


# TODO: Add a SpatiaLite geometry column
# eg. https://docs.datasette.io/en/stable/spatialite.html
#
# Spatial indexing latitude/longitude columns
# Here's a recipe for taking a table with existing latitude and longitude columns, adding a SpatiaLite POINT geometry column to that table, populating the new column and then populating a spatial index:
# 
# import sqlite3
# conn = sqlite3.connect('museums.db')
# # Lead the spatialite extension:
# conn.enable_load_extension(True)
# conn.load_extension('/usr/local/lib/mod_spatialite.dylib')
# # Initialize spatial metadata for this database:
# conn.execute('select InitSpatialMetadata(1)')
# # Add a geometry column called point_geom to our museums table:
# conn.execute("SELECT AddGeometryColumn('museums', 'point_geom', 4326, 'POINT', 2);")
# # Now update that geometry column with the lat/lon points
# conn.execute('''
#     UPDATE museums SET
#     point_geom = GeomFromText('POINT('||"longitude"||' '||"latitude"||')',4326);
# ''')
# # Now add a spatial index to that column
# conn.execute('select CreateSpatialIndex("museums", "point_geom");')
# # If you don't commit your changes will not be persisted:
# conn.commit()
# conn.close()
