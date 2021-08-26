import argparse
from pathlib import Path

from pyproj import Geod
from database import db
from shapely.prepared import prep
from shapely.wkb import load
import shapely.geometry as geometry

# specify a named ellipsoid
geod = Geod(ellps="WGS84")

parser = argparse.ArgumentParser(
    description="Compare the area of areas given as polygons with the total area of the corresponding voronoi polygons in the database"
)
parser.add_argument("wkb", nargs="+", type=Path)
args = parser.parse_args()

# args.wkb = [
#     Path("../plotting/{}.wkb".format(t))
#     for t in {
#         "Alaska",
#         "Haida Nation islands",
#         "California",
#         "Louisiana",
#         "Florida",
#         "Newfoundland",
#         "Baja California Sur",
#         "Cuba",
#         "Amazonas",
#         "Paraguay",
#         "Rio Negro (ARG)",
#         "Tierra del Fuego (Isla Grande)",
#     }
# ]

shapes = {}
for file in args.wkb:
    shapes[file.stem] = load(file.open("rb"))

DATABASE, TABLES = db()

print("Name", "Geodesic Area", "Voronoi Area", "Population Count", sep="\t")
for name, shape in shapes.items():
    east, south, west, north = shape.bounds

    area, _ = geod.geometry_area_perimeter(shape)

    candidate_points = DATABASE.execute(
        f"""SELECT node_id, min(latitude), min(longitude), max(popdensity), sum(area), sum(population_capacity)
        FROM nodes LEFT JOIN ecology ON node_id = node
        WHERE latitude < {north:} AND {south:} < latitude AND {east:} < longitude AND longitude < {west:} AND node_id > 100000000
        GROUP BY node_id"""
    ).all()

    shape = prep(shape)
    db_area = 0.0
    db_population = 0.0
    for (
        node_id,
        latitude,
        longitude,
        popdensity,
        voronoi_area,
        popcap,
    ) in candidate_points:
        if shape.contains(geometry.Point(longitude, latitude)):
            db_area += voronoi_area or 0.0
            db_population += popcap or 0.0

    print(name, area / 1000000, db_area, db_population, sep="\t")
