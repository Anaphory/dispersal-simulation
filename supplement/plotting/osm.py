#!/usr/bin/python2
# -*- encoding: utf-8 -*-
from pathlib import Path

import overpy
from shapely.prepared import prep
from shapely.wkb import load, dump
import shapely.geometry as geometry
from shapely.ops import linemerge, unary_union, polygonize
import cartopy.io.shapereader as shpreader

api = overpy.Overpass()

areas = {
    "Alaska": 1116270,
    "Haida Nation islands": 9585706,
    "California": 165475,
    "Louisiana": 224922,
    "Florida": 162050,
    "Newfoundland": 9467325,
    "Baja California Sur": 2589611,
    "Cuba": 307833,
    "Amazonas": 332476,
    "Paraguay": 287077,
    "Rio Negro (ARG)": 153547,
    "Tierra del Fuego (Isla Grande)": 163703,
}
shapes = {}

for area, id in areas.items():
    # Load and cache all prerequested areas
    fname = (Path(__file__).parent / "{:s}.wkb".format(area))
    try:
        shapes[area] = load(fname.open("rb"))
    except FileNotFoundError:
        print(u"Area {:s} not found in cache.".format(area))
        query = """[out:json][timeout:25];
        rel({:d});
        out body;
        >;
        out skel qt; """.format(
            id
        )
        result = api.query(query)

        # Convert ways to linstrings
        lss = []

        for ii_w, way in enumerate(result.ways):
            ls_coords = []

            for node in way.nodes:
                # create a list of node coordinates
                ls_coords.append((node.lon, node.lat))

            # create a LineString from coords
            lss.append(geometry.LineString(ls_coords))

        merged = linemerge(lss)  # merge LineStrings
        borders = unary_union(merged)  # linestrings to a MultiLineString
        polygons = list(polygonize(borders))
        shapes[area] = geometry.MultiPolygon(polygons)
        dump(shapes[area], fname.open("wb"))
        # FIXME: What is the canonical location for cached program supplement
        # data? Should I download it at install-time, when I'm copied somewhere
        # with write access that might never come back?


def contains(shape, coordinates):
    return shape.contains(geometry.Point(*coordinates))


land_shp_fname = shpreader.natural_earth(
    resolution="50m", category="physical", name="land"
)

land_geom = unary_union(
    [
        record.geometry
        for record in shpreader.Reader(land_shp_fname).records()
        if record.attributes.get("featurecla") != "Null island"
    ]
)

LAND = prep(land_geom)


def is_land(xy):
    return LAND.contains(geometry.Point(*xy))
