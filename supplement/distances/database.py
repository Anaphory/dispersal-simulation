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
        sqlalchemy.Column("travel_time", sqlalchemy.Float),
        sqlalchemy.Column("flat_distance", sqlalchemy.Float),
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
