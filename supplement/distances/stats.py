import numpy

from sqlalchemy import func
from sqlalchemy import select

from database import db

DATABASE, TABLES = db()

(n,) = DATABASE.execute(
    select(func.count(TABLES["nodes"].c.node_id)).where(
        TABLES["nodes"].c.node_id < 100000000
    )
).fetchone()
print(r"\newcommand{\rivernodes}{", n, "}")

(n,) = DATABASE.execute(
    select(func.count(TABLES["nodes"].c.node_id)).where(
        TABLES["nodes"].c.node_id > 100000000
    )
).fetchone()
print(r"\newcommand{\territorynodes}{", n, "}")

dists = numpy.array(
    DATABASE.execute(
        select(TABLES["edges"].c.travel_time / 3600.).where(
            (TABLES["edges"].c.node1 != TABLES["edges"].c.node2)
            & (TABLES["edges"].c.source == "sea")
        )
    ).fetchall()
)
print(
    r"\newcommand{\seadists}{",
    numpy.quantile(dists, [0.05, 0.5, 0.95]),
    numpy.mean(dists),
    "}",
)

dists = numpy.array(
    DATABASE.execute(
        select(TABLES["edges"].c.travel_time / 3600.).where(
            (TABLES["edges"].c.node1 != TABLES["edges"].c.node2)
            & (TABLES["edges"].c.travel_time < numpy.inf)
            & (TABLES["edges"].c.source == "grid")
        )
    ).fetchall()
)
print(
    r"\newcommand{\griddists}{",
    numpy.quantile(dists, [0.05, 0.5, 0.95]),
    numpy.mean(dists),
    "}",
)

dists = numpy.array(
    DATABASE.execute(
        select(TABLES["edges"].c.travel_time / 3600.).where(
            (TABLES["edges"].c.node1 > 100000000)
            & (TABLES["edges"].c.node2 > 100000000)
            & (TABLES["edges"].c.node1 != TABLES["edges"].c.node2)
            & (TABLES["edges"].c.travel_time < numpy.inf)
            & (TABLES["edges"].c.source == "grid")
        )
    ).fetchall()
)
print(
    r"\newcommand{\hexdists}{",
    numpy.quantile(dists, [0.05, 0.5, 0.95]),
    numpy.mean(dists),
    "}",
)

dists = numpy.array(
    DATABASE.execute(
        select(TABLES["edges"].c.travel_time / 3600.).where(
            (TABLES["edges"].c.node1 != TABLES["edges"].c.node2)
            & (TABLES["edges"].c.source == "river")
        )
    ).fetchall()
)
print(
    r"\newcommand{\riverdists}{",
    numpy.quantile(dists, [0.05, 0.5, 0.95]),
    numpy.mean(dists),
    "}",
)
