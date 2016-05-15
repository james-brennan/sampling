#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""

  Main routine

"""
from orbit_calculation import *
from swath_geom import *

import ephem
import datetime
import numpy as np
#import matplotlib.pyplot as plt

def main(name, line1, line2, orbital_filename, fov, swath_width):
    """
    Steps
    1. generate shapefile to save orbits
    2. Run generator which returns per orbital swath
        1. Create polygon per orbital swath
        2. Append to shapefile
    3. Save shapefile
    """
    satellite = ephem.readtle(name, line1, line2)
    

    # Landsat 8
    #name = "Landsat8"
    #line1="1 39084U 13008A   16051.82349873  .00000188  00000-0  51829-4 0  9999"
    #line2="2 39084  98.1988 123.2603 0001265  89.4360 270.6984 14.57110027160810"
    #LD8 = ephem.readtle(name, line1, line2)
    
    sun = ephem.Sun()
    fov = np.radians(fov) # convert to radians

    """
    Make pandas dataframe to store swath information
    """
    import pandas as pd
    data = {"DateTime": [],
            "orbit_id":[], "ground_lat": [], 
            "ground_lon": [], "swath_width": []}
    swaths = pd.DataFrame(data)
    swaths.set_index(keys="DateTime")
    # generate shapefile

    orbit_id = 0
    # need to do splitted by hemisphere unfortunately..
    for orbit in make_an_orbit(satellite):
        #import pdb; pdb.set_trace()
        if len(orbit) > 1:
            """
            So worth doing processing on orbit...

            """
            print(orbit[0].datetime)
            for overpass in orbit:
                overpass.only_daytime_overpasses(sun)
                overpass.derive_swath_width(fov)
            """
            Create a tempoary dataframe for this orbit
            """
            """
            Fix swath widths for now
            """
            #t = [o.swath_width = swath_width for o in orbit]
            epoch = datetime.datetime(1970, 1, 1)
            #import pdb; pdb.set_trace()
            tmp_d = {"DateTime": [(o.datetime - epoch).total_seconds() for o in orbit],
                     "orbit_id": orbit_id * np.ones(len(orbit)),
                     "ground_lat": [o.lat for o in orbit],
                     "ground_lon": [o.long for o in orbit],
                     "swath_width":  swath_width * np.ones(len(orbit))}
            tmp = pd.DataFrame(tmp_d)
            tmp.set_index(keys="DateTime")
            #import pdb; pdb.set_trace()
            orbit_id +=1 
            """
            Append to main dataframe
            """
            swaths = swaths.append(tmp)
            #swaths.set_index(keys="DateTime")

    """
    Save the DataFrame to a file
    """
    swaths = swaths.set_index(keys="DateTime")
    #swaths.set_index(keys="DateTime")
    #import pdb; pdb.set_trace()
    swaths.to_csv(orbital_filename, header=False)


if __name__ == "__main__":

    name = "TERRA"
    line1 = "1 25994U 99068A   16048.43680378  .00000258  00000-0  67198-4 0  9999"
    line2 = "2 25994  98.1982 124.4247 0001352 105.3907 254.7441 14.57126067859938"
    filename = "TERRA_orbits.csv"
    fov = 110.0# in degrees
    swath_width = 2330 * 1e3
    main(name, line1, line2, filename, fov, swath_width)

    name = "AQUA"
    line1 = "1 27424U 02022A   16057.87908387  .00000304  00000-0  77681-4 0  9994"
    line2 = "2 27424  98.2234   0.0918 0001365  95.1594 344.4495 14.57097852734910"
    filename = "AQUA_orbits.csv"
    fov = 110.0# in degrees
    swath_width = 2330 * 1e3
    main(name, line1, line2, filename, fov, swath_width)

    name = "Landsat8_OLI"
    line1="1 39084U 13008A   16051.82349873  .00000188  00000-0  51829-4 0  9999"
    line2="2 39084  98.1988 123.2603 0001265  89.4360 270.6984 14.57110027160810"
    filename = "Landsat8_OLI_orbits.csv"
    fov = 15# in degrees
    swath_width = 185* 1e3
    main(name, line1, line2, filename, fov, swath_width)

    name = "MERIS"
    line1="1 27386U 02009A   16058.06144687  .00000037  00000-0  26298-4 0  9994"
    line2="2 27386  98.2990 115.3795 0000852  79.1696 280.9586 14.37857112732659"
    filename = "MERIS_orbits.csv"
    fov = 68.5# in degrees
    swath_width = 1150 * 1e3
    main(name, line1, line2, filename, fov, swath_width)

    name = "SENTINEL-3A-OLCI"
    line1="1 41335U 16011A   16132.45512502 -.00000044  00000-0  00000+0 0  9990"
    line2="2 41335  98.6273 199.6294 0001406 112.0638 248.0700 14.26735166 12077"
    filename = "SENTINEL-3A-OLCI_orbits.csv"
    fov = 68.5# in degrees
    swath_width = 1270  * 1e3
    main(name, line1, line2, filename, fov, swath_width)

    name = "SENTINEL-2A_MSI"
    line1="1 40697U 15028A   16132.43135143  .00000139  00000-0  69694-4 0  9990"
    line2="2 40697  98.5685 207.1695 0001293  97.8190 262.3139 14.30816868 46245"
    filename = "SENTINEL-2A_MSI_orbits.csv"
    fov = 68.5# in degrees
    swath_width = 290  * 1e3
    main(name, line1, line2, filename, fov, swath_width)
