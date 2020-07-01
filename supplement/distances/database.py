import typing as t

import sqlalchemy
from sqlalchemy import event
from sqlalchemy.engine import Engine

from h3 import h3

@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection: sqlalchemy.engine.Connectable, connection_record: t.Any) -> None:
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()


def db(file: str = "sqlite:///hexbins.sqlite") -> t.Tuple[sqlalchemy.engine.Connectable, t.Dict[str, sqlalchemy.Table]]:
    class H3Index(sqlalchemy.types.TypeDecorator): # type: ignore
        impl = sqlalchemy.Integer

        def process_bind_param(self, value: str, dialect: t.Any) -> sqlalchemy.Integer:
            return h3.string_to_h3(value)

        def process_result_value(self, value: sqlalchemy.Integer, dialect: t.Any) -> str:
            return h3.h3_to_string(value)

    engine = sqlalchemy.create_engine(file)
    metadata = sqlalchemy.MetaData()
    hex = sqlalchemy.Table(
        'hex', metadata,
        sqlalchemy.Column('hexbin', H3Index, primary_key=True),
        sqlalchemy.Column('longitude', sqlalchemy.Float),
        sqlalchemy.Column('latitude', sqlalchemy.Float),
        sqlalchemy.Column('vlongitude', sqlalchemy.Float),
        sqlalchemy.Column('vlatitude', sqlalchemy.Float),
        sqlalchemy.Column('habitable', sqlalchemy.Boolean()),
    )
    dist = sqlalchemy.Table(
        'dist', metadata,
        sqlalchemy.Column(
            'hexbin1', H3Index,
            sqlalchemy.ForeignKey(hex.c.hexbin),
            primary_key=True),
        sqlalchemy.Column(
            'hexbin2', H3Index,
            sqlalchemy.ForeignKey(hex.c.hexbin),
            primary_key=True),
        sqlalchemy.Column('distance', sqlalchemy.Float),
        sqlalchemy.Column('flat_distance', sqlalchemy.Float),
        sqlalchemy.Column('source', sqlalchemy.Integer),
    )
    eco = sqlalchemy.Table(
        'eco', metadata,
        sqlalchemy.Column(
            'hexbin', H3Index,
            sqlalchemy.ForeignKey(hex.c.hexbin), primary_key=True),
        sqlalchemy.Column('ecoregion', sqlalchemy.Integer(), primary_key=True),
        sqlalchemy.Column('frequency', sqlalchemy.Float),
    )
    reach = sqlalchemy.Table(
        'reach', metadata,
        sqlalchemy.Column('reach_id', sqlalchemy.Integer(), primary_key=True),
        sqlalchemy.Column('record_index', sqlalchemy.Integer()),
        sqlalchemy.Column('next_down', sqlalchemy.Integer(), sqlalchemy.ForeignKey('reach.reach_id')),
        sqlalchemy.Column("travel_downstream", sqlalchemy.Float()),
        sqlalchemy.Column("travel_upstream", sqlalchemy.Float()),
    )
    flows = sqlalchemy.Table(
        'flows', metadata,
        sqlalchemy.Column(
            'hexbin', H3Index,
            sqlalchemy.ForeignKey(hex.c.hexbin), primary_key=True),
        sqlalchemy.Column(
            'reach_id', sqlalchemy.Integer(),
            sqlalchemy.ForeignKey(reach.c.reach_id), primary_key=True),
        sqlalchemy.Column('start', sqlalchemy.Float()),
        sqlalchemy.Column('end', sqlalchemy.Float()),
    )
    metadata.create_all(engine)
    return engine, {
        "hex": hex, "dist": dist, "eco": eco, "reach": reach, "flows": flows}
