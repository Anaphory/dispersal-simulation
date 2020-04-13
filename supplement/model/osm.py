import overpy
import shapely.geometry as geometry
from shapely.wkb import load, dump
from shapely.ops import linemerge, unary_union, polygonize

api = overpy.Overpass()

areas = {
    "Alaska": 1116270,
}
shapes = {}

for area, id in areas.items():
    try:
        shapes[area] = load(open("{:s}.wkb".format(area), "rb"))
    except FileNotFoundError:
        print("Area {:s} not found in cache.".format(area))
        query = """[out:json][timeout:25];
        rel({:d});
        out body;
        >;
        out skel qt; """.format(id)
        result = api.query(query)

        lss = [] #convert ways to linstrings

        for ii_w,way in enumerate(result.ways):
            ls_coords = []

            for node in way.nodes:
                ls_coords.append((node.lon,node.lat)) # create a list of node coordinates

            lss.append(geometry.LineString(ls_coords)) # create a LineString from coords


        merged = linemerge([*lss]) # merge LineStrings
        borders = unary_union(merged) # linestrings to a MultiLineString
        polygons = list(polygonize(borders))
        shapes[area] = geometry.MultiPolygon(polygons)
        dump(shapes[area], open("{:s}.wkb".format(area), "wb"))

def contains(shape, coordinates):
    return shape.contains(geometry.Point(*coordinates))


