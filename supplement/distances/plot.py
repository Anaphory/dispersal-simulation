import numpy
import sqlalchemy
from sqlalchemy.sql import func
import matplotlib.pyplot as plt
import sys

from database import db
from ecoregions import TC
from raster_data import ecoregion_tile_from_geocoordinates


def plot_distances(engine: sqlalchemy.engine.Connectable, dist: sqlalchemy.Table) -> None:
    distances = sqlalchemy.select([dist.c.flat_distance, dist.c.distance]).where(
        dist.c.source == 0)
    x, y = zip(*engine.execute(distances))
    plt.scatter(x, y, marker='x', s=40, alpha=0.2)


def plot_locations(engine: sqlalchemy.engine.Connectable, t_hex: sqlalchemy.Table) -> None:
    d = sqlalchemy.select(
        [t_hex.c.vlongitude, t_hex.c.vlatitude, t_hex.c.habitable]
    ).where(
        t_hex.c.vlatitude is not None
    )
    x, y, h = zip(*engine.execute(d))
    h = numpy.array(h, dtype=bool)
    plt.scatter(
        x, y, marker='o', alpha=0.5, c=numpy.array(list("rg"))[1*h])
    er = ecoregion_tile_from_geocoordinates(-165, 60).read(1)
    plt.imshow(TC[er], extent=(-180, -150, 50, 70))


def plot_sampled(engine: sqlalchemy.engine.Connectable, t_hex: sqlalchemy.Table) -> None:
    d = sqlalchemy.select(
        [t_hex.c.vlongitude, t_hex.c.vlatitude, func.min(t_dist.c.distance), t_hex.c.habitable]
    ).where(
        t_hex.c.hexbin == t_dist.c.hexbin1
    ).where(
        t_dist.c.distance > 0
    ).group_by(
        t_hex.c.hexbin
    )
    x, y, s, h = zip(*engine.execute(d))
    s = numpy.array(s)
    h = numpy.array(h)
    plt.scatter(
        x, y, marker='o', s=s/3600, alpha=0.3, c='r')
    er = ecoregion_tile_from_geocoordinates(-165, 60).read(1)
    plt.imshow(-TC[er], extent=(-180, -150, 50, 70))


if __name__ == '__main__':
    # FIXME: Use Argparser instead
    import sys
    engine, tables = db(sys.argv[1])
    t_hex = tables["hex"]
    t_dist = tables["dist"]
    t_eco = tables["eco"]
    plot_distances(engine, t_dist)
    plt.show()
    plot_sampled(engine, t_hex)
    plt.show()
    plot_locations(engine, t_hex)
    plt.show()
