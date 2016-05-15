
# run everything
#python main.py
#echo "Finished orbits calculation Now running counting"
#echo date
# run the counting
#./calculateOverpasses TERRA_orbits.csv terra_counts.tif
#./calculateOverpasses AQUA_orbits.csv aqua_counts.tif
#./calculateOverpasses Landsat8_OLI_orbits.csv landsat8_OLI_counts.tif
#./calculateOverpasses MERIS_orbits.csv meris_counts.tif
#./calculateOverpasses SENTINEL-3A-OLCI_orbits.csv sentinel3a_OLCI_counts.tif
#./calculateOverpasses SENTINEL-2A_MSI_orbits.csv sentinel2a_MSI_counts.tif

nohup ./samplingModel aqua_counts.tif aqua_days.tif &
nohup ./samplingModel terra_counts.tif terra_days.tif &
nohup ./samplingModel landsat8_OLI_counts.tif L8_OLI_days.tif &
nohup ./samplingModel meris_counts.tif meris_days.tif &
nohup ./samplingModel sentinel3a_OLCI_counts.tif sentinel3A_OLCI_days.tif &
nohup ./samplingModel sentinel2a_MSI_counts.tif sentinel2A_MSI_days.tif &
