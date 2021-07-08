import typing as t

import cartopy.geodesic as geodesic
import shapely.geometry as sgeom
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader
from shapely.prepared import prep

LonLat = t.Tuple[float, float]

# Define some constants
GEODESIC: geodesic.Geodesic = geodesic.Geodesic()

BBOX = sgeom.box(-168.956, -55.852, -17.556, 83.696)  # minx, miny, maxx, maxy

ALL_LAND = unary_union(
    [
        record.geometry
        for record in shpreader.Reader(
            shpreader.natural_earth(resolution="10m", category="physical", name="land")
        ).records()
        if record.attributes.get("featurecla") != "Null island"
    ]
).difference(
    unary_union(
        list(
            shpreader.Reader(
                shpreader.natural_earth(
                    resolution="10m", category="physical", name="lakes"
                )
            ).geometries()
        )
    )
)
LAND = BBOX.intersection(ALL_LAND)
PLAND = prep(LAND)
DEFINITELY_INLAND = prep(LAND.buffer(0.03))
