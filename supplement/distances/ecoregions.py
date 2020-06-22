import zipfile
import typing as t
from pathlib import Path

import shapefile
import shapely.geometry as sgeom
from shapely.prepared import prep


class Ecoregions:
    cache = None
    @classmethod
    def ecoregions(cls):
        """
        >>> eco = Ecoregions.ecoregions()
        >>> eco.numRecords
        847
        >>> eco.fields
        [('DeletionFlag', 'C', 1, 0), ['OBJECTID', 'N', 32, 10], ['ECO_NAME', 'C', 150, 0], ['BIOME_NUM', 'N', 32, 10], ['BIOME_NAME', 'C', 254, 0], ['REALM', 'C', 254, 0], ['ECO_BIOME_', 'C', 254, 0], ['NNH', 'N', 11, 0], ['ECO_ID', 'N', 11, 0], ['SHAPE_LENG', 'N', 32, 10], ['SHAPE_AREA', 'N', 32, 10], ['NNH_NAME', 'C', 64, 0], ['COLOR', 'C', 7, 0], ['COLOR_BIO', 'C', 7, 0], ['COLOR_NNH', 'C', 7, 0], ['LICENSE', 'C', 64, 0]]

        """
        if cls.cache is not None:
            return cls.cache
        zipshape = zipfile.ZipFile(
            (Path(__file__).parent /
             "../ecoregions/Ecoregions2017.zip").open("rb"))
        shape = shapefile.Reader(
            shp=zipshape.open("Ecoregions2017.shp"),
            shx=zipshape.open("Ecoregions2017.shx"),
            dbf=zipshape.open("Ecoregions2017.dbf"),
            encoding='latin-1'
        )
        cls.cache = shape
        return shape

    def __init__(self, mask: t.Optional[sgeom.Polygon] = None):
        """

        >>> eco = Ecoregions()
        >>> eco.records[508][1]
        'Tocantins/Pindare moist forests'

        """
        shp = self.ecoregions()
        self.n = shp.numRecords
        self.records = {int(shp.record(eco)[0]): shp.record(eco)
                        for eco in range(self.n)}
        self.mask = mask
        self._b = None

    @property
    def boundaries(self):
        if self._b is not None:
            return self._b

        self._b = []
        for b, boundary in enumerate(shp.shapes()):
            if self.mask is None:
                self._b.append(prep(sgeom.shape(boundary)))
            else:
                # 205 is the ice shield, which intersects everything.
                if b != 205 and sgeom.box(*boundary.bbox).intersects(self.mask):
                    s = sgeom.shape(boundary).buffer(0)
                    self._b.append(prep(s & self.mask))
                else:
                    # The bounding box of the feature is disjoint with the
                    # mask, nothing to do.
                    self._b.append(None)

        return self._b


    def at_point(self, point):
        """Get the ecoregion for the given location

        >>> ECOREGION = Ecoregions()
        >>> r = ECOREGION.at_point(sgeom.Point(-38.282, -8.334))
        >>> r
        84
        >>> ECOREGION.record_getter(r)
        Record #84: [86.0, 'Caatinga', 2.0, 'Tropical & Subtropical Dry Broadleaf Forests', 'Neotropic', 'NO02', 4, 525, 113.696454941, 60.155581182, 'Nature Imperiled', '#81B50A', '#CCCD65', '#EE1E23', 'CC-BY 4.0']
        """
        for b, boundary in enumerate(self.boundaries):
            if boundary and boundary.contains(point):
                return b
        return None

ECOREGIONS = Ecoregions()

