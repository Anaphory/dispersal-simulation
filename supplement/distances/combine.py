import json
import sqlite3
import collections
from pathlib import Path

import h3
import numpy
import networkx

import matplotlib.pyplot as plt

INF = float("inf")

try:
    Path("plot.sqlite").unlink()
except FileNotFoundError:
    pass

conn = sqlite3.connect("plot.sqlite")

conn.execute(
    """
    CREATE TABLE hex (hexbin INTEGER PRIMARY KEY, longitude REAL, latitude REAL)
"""
)
conn.execute(
    """
    CREATE TABLE dist (hexbin1 INTEGER, hexbin2 INTEGER, distance REAL, PRIMARY KEY (hexbin1, hexbin2))
"""
)
conn.execute(
    """
    CREATE TABLE eco (hexbin INTEGER, ecoregion INTEGER, area_m2 REAL, PRIMARY KEY (hexbin, ecoregion))
"""
)

g = networkx.Graph()

dist_table = {}

hexes = set()
for file in Path().glob("hexdist*.json"):
    with file.open() as dist_json:
        dist = json.load(dist_json)
        for start, connections in dist.items():
            start = int(start)
            hexes.add(start)
            for end, distance in connections.items():
                end = int(end)
                hexes.add(end)
                dist_table[start, end] = min(
                    dist_table.get((start, end), INF), distance
                )

conn.executemany(
    """
    INSERT INTO dist VALUES (?, ?, ?)
""",
    ((h1, h2, d) for (h1, h2), d in dist_table.items()),
)


eco_table = collections.defaultdict(collections.Counter)
for file in Path().glob("hex_to_eco*.json"):
    with file.open() as eco_json:
        eco = json.load(eco_json)
        for start, eco in eco.items():
            start = int(start)
            hexes.add(start)
            for type, amount in eco.items():
                type = int(type)
                eco_table[start][type] += amount

conn.executemany(
    """
    INSERT INTO eco VALUES (?, ?, ?)
""",
    ((h, t, a) for h, regions in eco_table.items() for t, a in regions.items()),
)

hex_table = {}
for hex in hexes:
    lat, lon = h3.api.basic_int.h3_to_geo(hex)
    hex_table[hex] = lon, lat

conn.executemany(
    """
    INSERT INTO hex VALUES (?, ?, ?)
""",
    ((h, x, y) for h, (x, y) in hex_table.items()),
)

conn.commit()

g.add_weighted_edges_from((h1, h2, d) for (h1, h2), d in dist_table.items())
networkx.draw_networkx_edges(
    g,
    hex_table,
    edge_color=1
    / numpy.array(list(networkx.get_edge_attributes(g, "weight").values())),
    edge_cmap=plt.cm.Blues,
)
networkx.draw_networkx_nodes(
    g,
    hex_table,
    node_size=1,
)
plt.show()
