for zip in ../elevation/GMTED2010/*.zip
do
    file=$(unzip -l $zip | grep -o '[^ ]*med150.tif$')
    echo $file
    extent=$(gdalinfo /vsizip/$zip/$file | grep 'Lower Left\|Upper Right' | sed -E 's/[^\(]+\(([^,]+),([^\)]+)\).+/\1 \2/' | paste --serial)
    gdalsrsinfo -o wkt /vsizip/$zip/$file > GMTED.wkt
    echo gdal_rasterize -l Ecoregions2017 -ts 7200 4800 -a OBJECTID -a_nodata 999 -a_srs GMTED.wkt -ot UInt16 -of GTiff /vsizip/Ecoregions2017.zip -te $extent /tmp/ECOREGIONS-${file}
    gdal_rasterize -l Ecoregions2017 -ts 7200 4800 -a OBJECTID -a_nodata 999 -a_srs GMTED.wkt -ot UInt16 -of GTiff /vsizip/Ecoregions2017.zip -te $extent ECOREGIONS-${file}
done
