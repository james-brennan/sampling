#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""

  Swath_geom.py


  Does some of the work of generating swath shapefiles from
  a ground track

  Basically take all of the overpasses as turn into a polygon

  --- current issue above about 70°
  =================================
        this shouldn't matter as we're interested in BA

"""
from osgeo import ogr
from osgeo import osr
import numpy as np
from shapely.geometry.polygon import LinearRing
from osgeo import ogr
from shapely.geometry import Polygon






def make_swth_polygon_shift(overpasses):
    """

    Sort of works but issue about wrapping longitude atm

    solution:
        ONLY PROVIDE dependent on meridian eg pos or neg longitude
        need a function which will split these...

        Also only provide seperately for ascending and descending
    """

    """
    Limit to within latitudes...
    70s to 70N
    """
    #import pdb; pdb.set_trace()
    limited_pass = [o for o in overpasses if ((o.lat > -70) & (o.lat < 70) )]
    """
    Split coords into leftmost and rightmost

    think this needs to be cleverer...
        eg this will change if satelitte is ascending or descending

    idea 1: check bottom and top for both boundary longs
        the left one should be less for both

    """
    import pdb; pdb.set_trace()
    first_coords = np.array([[o.lats[0], o.lons[0]] for o in limited_pass])
    second_coords = np.array([[o.lats[1], o.lons[1]] for o in limited_pass])
    from shapely.geometry import LineString
    all_coords = np.concatenate((first_coords, second_coords[::-1]))
    # make a linestring from these coords
    linestring = LineString(all_coords)
    # shift longitude
    import shift_geom
    l2 = shift_geom.shift(linestring)
    x = all_coords[:,1]
    y = all_coords[:,0]
    xy = np.dstack((x,y))
    #import pdb; pdb.set_trace()
    #LR = LinearRing(xy[0])
    LR = LineString(xy[0])
    LR2 = shift_geom.shift(LR)
    poly = Polygon(LR2)
    return LR2








def make_two_swath_linestrings(overpasses):
    """

    Sort of works but issue about wrapping longitude atm

    solution:
        ONLY PROVIDE dependent on meridian eg pos or neg longitude
        need a function which will split these...

        Also only provide seperately for ascending and descending
    """

    """
    Limit to within latitudes...
    70s to 70N
    """
    #import pdb; pdb.set_trace()
    limited_pass = [o for o in overpasses if ((o.lat > -70) & (o.lat < 70) )]
    """
    Split coords into leftmost and rightmost

    think this needs to be cleverer...
        eg this will change if satelitte is ascending or descending

    idea 1: check bottom and top for both boundary longs
        the left one should be less for both

    """
    #import pdb; pdb.set_trace()
    first_coords = np.array([[o.lats[0], o.lons[0]] for o in limited_pass])
    second_coords = np.array([[o.lats[1], o.lons[1]] for o in limited_pass])
    from shapely.geometry import LineString, MultiLineString
    from shapely.ops import cascaded_union
    x1 = first_coords[:,1]
    y1 = first_coords[:,0]
    xy1 = np.dstack((x1,y1))
    #import pdb; pdb.set_trace()
    #LR = LinearRing(xy[0])
    line1 = LineString(xy1[0])
    x1 = second_coords[:,1]
    y1 = second_coords[:,0]
    xy2 = np.dstack((x1,y1))
    line2 = LineString(xy2[0])

    mls = MultiLineString([line1, line2])

    return mls



#ogr2ogr -f "ESRI Shapefile" -overwrite test2.shp test.shp -sql "SELECT * FROM testing" -skipfailures -nlt MULTIPOLYGON


def make_swth_polygon2(overpasses):
    """

    Sort of works but issue about wrapping longitude atm

    solution:
        ONLY PROVIDE dependent on meridian eg pos or neg longitude
        need a function which will split these...

        Also only provide seperately for ascending and descending
    """

    """
    Limit to within latitudes...
    70s to 70N
    """
    #import pdb; pdb.set_trace()
    limited_pass = [o for o in overpasses if ((o.lat > -70) & (o.lat < 70) )]
    """
    Split coords into leftmost and rightmost

    think this needs to be cleverer...
        eg this will change if satelitte is ascending or descending

    idea 1: check bottom and top for both boundary longs
        the left one should be less for both

    """
    #import pdb; pdb.set_trace()
    first_coords = np.array([[o.lats[0], o.lons[0]] for o in limited_pass])
    second_coords = np.array([[o.lats[1], o.lons[1]] for o in limited_pass])

    west1 =[]
    east1 = []
    west2 = []
    east2 = []
    if np.any(first_coords[:, 1] > 0):
        # some are in the east and some in the west...
        # bollocks
        # need to make two seperate polygons...
        east_idx1 = first_coords[:, 1] > 0
        west_idx1 = first_coords[:, 1] < 0
    if np.any(second_coords[:, 1] > 0):
        # some are in the east and some in the west...
        # bollocks
        # need to make two seperate polygons...
        east_idx2 = second_coords[:, 1] > 0
        west_idx2 = second_coords[:, 1] < 0

    # now need to do?
    import pdb; pdb.set_trace()


    """
    Or anywhere it is negative reflect to make positive?
    """
    #west = first_coords[:, 1] < 0
    #first_coords[west, 1] = abs(first_coords[west, 1])
    #west = second_coords[:, 1] < 0
    #second_coords[west, 1] = abs(second_coords[west, 1])
    #first_coords[:, 1] = (((first_coords[:, 1] - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
    #first_coords[:, 1] += 180
    #second_coords[:, 1] = (((second_coords[:, 1] - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin

    """
    How to solve negative wrap around
    1. Move westly points (negative long) + 360
    2. Move whats not moved (pos long) + 360

    """


    key = {0: first_coords, 1: second_coords}
    #first_coords[:, 1] = (first_coords[:, 1] + 180) % 360
    #second_coords[:, 1] = (second_coords[:, 1] + 180) % 360
    top_longs = (first_coords[0][1], second_coords[0][1]   )
    bottom_longs =  (first_coords[-1][1], second_coords[-1][1] )
    #import pdb; pdb.set_trace()
    """
    Determine which of the set is nearest...
    """
    assert np.argmin(top_longs) == np.argmin(bottom_longs) # should be the same
    nearest_to_left = np.argmin(top_longs)
    nearest_to_right = np.argmax(top_longs)


    left_coords = key[nearest_to_left]
    right_coords = key[nearest_to_right]
    # need to reverse so that they draw in right order!
    all_coords = np.concatenate((left_coords, right_coords[::-1]))
    """
    LR needs to be specified in right order


    eg [bl to br, all up right, tr to tl, all down left, attach left[-1] to bl]

    """
    #import pdb; pdb.set_trace()
    """
    First two points are the bottom of left and right
    """
    #bl = (left_coords[0])
    #br = (right_coors[0])
    """
    last two are tops of both
    """
    #tl = (left_coords[-1])
    #tr = (right_coors[-1])
    #bottom_connect = [bl, br]
    #right_ascending = [rc for  rc in right_coors[1:-1]]
    #top_connect = [tr, tl]
    #left_descending = [lc for lc in reversed(left_coords[1:-1])]
    #connectall = [left_descending[-1], bl]
    # flatten points to a list
    #points = np.vstack([bottom_connect, right_ascending, top_connect, left_descending])
    """

    Shapely uses x,y so need long, lat

    """
    x = all_coords[:,1]
    y = all_coords[:,0]
    xy = np.dstack((x,y))
    #import pdb; pdb.set_trace()
    #LR = LinearRing(xy[0])
    poly = Polygon(xy[0])
    return poly










def make_swth_polygon(overpasses):
    """

    Sort of works but issue about wrapping longitude atm

    solution:
        ONLY PROVIDE dependent on meridian eg pos or neg longitude
        need a function which will split these...

        Also only provide seperately for ascending and descending
    """

    """
    Limit to within latitudes...
    70s to 70N
    """
    #import pdb; pdb.set_trace()
    limited_pass = [o for o in overpasses if ((o.lat > -70) & (o.lat < 70) )]
    """
    Split coords into leftmost and rightmost

    think this needs to be cleverer...
        eg this will change if satelitte is ascending or descending

    idea 1: check bottom and top for both boundary longs
        the left one should be less for both

    """
    #import pdb; pdb.set_trace()
    first_coords = np.array([[o.lats[0], o.lons[0]] for o in limited_pass])
    second_coords = np.array([[o.lats[1], o.lons[1]] for o in limited_pass])
    """
    Want to rescale longitude from -180..180 to 0..360
    NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin

    """
    import pylab as plt
    NewMax = 360
    NewMin = 0
    OldMin = -180
    OldMax = 180

    """
    Or anywhere it is negative reflect to make positive?
    """
    #west = first_coords[:, 1] < 0
    #first_coords[west, 1] = abs(first_coords[west, 1])
    #west = second_coords[:, 1] < 0
    #second_coords[west, 1] = abs(second_coords[west, 1])
    #first_coords[:, 1] = (((first_coords[:, 1] - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
    #first_coords[:, 1] += 180
    #second_coords[:, 1] = (((second_coords[:, 1] - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin

    """
    How to solve negative wrap around
    1. Move westly points (negative long) + 360
    2. Move whats not moved (pos long) + 360

    """


    key = {0: first_coords, 1: second_coords}
    #first_coords[:, 1] = (first_coords[:, 1] + 180) % 360
    #second_coords[:, 1] = (second_coords[:, 1] + 180) % 360
    top_longs = (first_coords[0][1], second_coords[0][1]   )
    bottom_longs =  (first_coords[-1][1], second_coords[-1][1] )
    #import pdb; pdb.set_trace()
    """
    Determine which of the set is nearest...
    """
    #assert np.argmin(top_longs) == np.argmin(bottom_longs) # should be the same
    nearest_to_left = np.argmin(top_longs)
    nearest_to_right = np.argmax(top_longs)


    left_coords = key[nearest_to_left]
    right_coords = key[nearest_to_right]
    # need to reverse so that they draw in right order!
    all_coords = np.concatenate((left_coords, right_coords[::-1]))
    """
    LR needs to be specified in right order


    eg [bl to br, all up right, tr to tl, all down left, attach left[-1] to bl]

    """
    #import pdb; pdb.set_trace()
    """
    First two points are the bottom of left and right
    """
    #bl = (left_coords[0])
    #br = (right_coors[0])
    """
    last two are tops of both
    """
    #tl = (left_coords[-1])
    #tr = (right_coors[-1])
    #bottom_connect = [bl, br]
    #right_ascending = [rc for  rc in right_coors[1:-1]]
    #top_connect = [tr, tl]
    #left_descending = [lc for lc in reversed(left_coords[1:-1])]
    #connectall = [left_descending[-1], bl]
    # flatten points to a list
    #points = np.vstack([bottom_connect, right_ascending, top_connect, left_descending])
    """

    Shapely uses x,y so need long, lat

    """
    x = all_coords[:,1]
    y = all_coords[:,0]
    xy = np.dstack((x,y))
    #import pdb; pdb.set_trace()
    #LR = LinearRing(xy[0])
    poly = Polygon(xy[0])
    return poly
