# reproject monthlies to lower resolution
res=2
while read file
do 
gdal_translate -tr $res $res $file ${res}_${file}
done < monthlys
# and do the lw mask
gdal_translate -tr $res $res cv_01_05_1km_uint16.tif ${res}_LW_mask.tif
