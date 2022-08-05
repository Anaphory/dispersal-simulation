for file in ../elevation/GMTED2010/*.tif
do
    echo $file
    basename=$(basename $file)
    extent=$(gdalinfo $file | grep 'Lower Left\|Upper Right' | sed -E 's/[^\(]+\(([^,]+),([^\)]+)\).+/\1 \2/' | paste --serial)
    gdalsrsinfo -o wkt $file > GMTED.wkt
    if [ -f ECOREGIONS-${basename} ]
    then
        echo ECOREGIONS-${basename}
    else
        echo gdal_rasterize -l Ecoregions2017 -ts 7200 4800 -a OBJECTID -a_nodata 999 -a_srs GMTED.wkt -ot UInt16 -of GTiff /vsizip/Ecoregions2017.zip -te $extent ECOREGIONS-${basename}
             gdal_rasterize -l Ecoregions2017 -ts 7200 4800 -a OBJECTID -a_nodata 999 -a_srs GMTED.wkt -ot UInt16 -of GTiff /vsizip/Ecoregions2017.zip -te $extent ECOREGIONS-${basename} &
    fi
done
