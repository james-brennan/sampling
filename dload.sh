
for x in $(seq -f "%02g" 01 12)
do
wget http://data.earthenv.org/cloud/MODCF_monthlymean_$x.tif
done
