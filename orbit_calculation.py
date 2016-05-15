#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""

Minimal working example...

"""
import utm
import ephem
import datetime
import numpy as np
import geopy
from geopy.distance import VincentyDistance



"""
Utility functions
"""
# from http://stackoverflow.com/questions/10688006/generate-a-list-of-datetimes-between-an-interval-in-python
def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

#class orbit(object):





"""

Method to make orbit
split when either descending or ascending and return


"""
def make_an_orbit(sensor_pyephem, get_west = False):
    overpasses = []
    eclipsed = False

    start = datetime.datetime(2016, 1, 1, 0,0,0,0)
    end =  datetime.datetime(2017, 1, 1, 0,0,0,0)
    delta = datetime.timedelta(seconds=10)
    """
    satelitte viewing day...
    """
    for time in perdelta(start, end, delta):
        ov = get_overpass(sensor_pyephem, time)
        overpasses.append(ov)
        """
        Theory:
        if the satelitte has been eclipsed it has travelled behind
        the earth and therefore over the top... or bottom
        """
        #if (ov.west == get_west) or (ov.eclipsed) or not ((ov.lat > -70) & (ov.lat < 70)):
        if (ov.eclipsed)  or  not ((ov.lat > -85) & (ov.lat < 80)):
            """
            Is shadowed by the earth...
            so good point to split and yield some overpasses
            """
            eclipsed = True
            # remove it
            overpasses.pop()
            yield overpasses
            #empty overpasses
            overpasses = []
            # change start time...
            start = time


"""
    Orbital overpass instance for a location, time etc
"""
class overpass(object):
    """
    Class to store overpass for a location
    """
    def __init__(self, lat, long, datetime, elevation, eclipsed, west):
        self.lat = lat
        self.long = long
        self.datetime = datetime
        self.daytime = None
        self.solarZenith = None
        self.sat_elevation = elevation # elevation of satellite
        self.eclipsed = eclipsed # if eclisped by the earth --> same as night?
        self.swath_width = None
        self.lats = None
        self.lons = None
        self.heading = None
        self.true_heading = None
        self.west = west

    def derive_swath_width(self, fov):
        """

        """
        R_e = 6378137.0 # earth radius [m]
        s = 0.5 * fov # half of field of view
        h = self.sat_elevation
        sin_alpha = np.arcsin( ( np.sin(s) * (R_e + h) ) / ( R_e )   ) - s
        alpha = np.sin(sin_alpha)
        swath = 2*((alpha/2*np.pi) * R_e) #metres
        self.swath_width =  swath

    def calculate_heading(self, dOverpass):
        """
        Not sure about an pyephem api method so
        just calculating heading by doing a
        difference over lat_0 -> lat_1 and lon_0 _> lat_1
        NOTE:
            delta_lat / delta_lon need to be small!!!
        """
        dlat = dOverpass.lat - self.lat
        dlon = dOverpass.long - self.long
        #assert dlat < 0.05 # make sure its small for now!
        #assert dlon < 0.05 # make sure its small for now!
        # incase provided backward
        # the direction is as vector from lat + dlat
        vlat = self.lat + dlat
        vlon = self.long + dlon
        """
        calculate heading!
        -- seems geopy does'nt do this yet

        https://gist.github.com/jeromer/2005586 has an answer which
        is re-used here...


        The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))

        """
        diffLong = np.radians(dOverpass.long  - self.long)
        lat1 = self.lat
        lat2 = dOverpass.lat
        x = np.sin(diffLong) * np.cos(lat2)
        y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1)
            * np.cos(lat2) * np.cos(diffLong))
        initial_bearing = np.arctan2(x, y)
        # make into a compass bearing
        initial_bearing = np.degrees(initial_bearing)
        self.true_heading  = initial_bearing
        compass_bearing = (initial_bearing + 360) % 360
        self.heading = compass_bearing


    def calculate_swath_edges(self):
        """
        returns min_lat, min_lon and max_lat,_max edges of swath
        actually seems easier to use geopy...
        """
        # use swath centre -- ground track point
        # 1. calculate bearing somehow... perpendicular to heading ...
        """

        Calculate two perpendiculars to the heading of the satellite

        perp1 = heading - 90
        perp1 = heading + 90
        # normalise to compass bearings!
        """
        perp1 = self.heading - 90
        perp2 = self.heading + 90
        perp1 = (perp1 + 360) % 360
        perp2 = (perp2 + 360) % 360
        """
        Do perp1
        """
        distance = self.swath_width / 2.0
        origin = geopy.Point(self.lat, self.long)
        destination1 = VincentyDistance(meters=distance).destination(origin, perp1)
        """
        Do perp2
        """
        destination2 = VincentyDistance(meters=distance).destination(origin, perp2)
        """ 
        Save as utm easting and northings...

        """
        #destination1 = utm.from_latlon(*destination1)
        #destination2 = utm.from_latlon(*destination2)
        self.lats = [destination1[0], destination2[0]]
        self.lons = [destination1[1], destination2[1]]         
        #import pdb; pdb.set_trace()
        #if destination1[1] < destination2[1]:
        #    # use as edges
        #    self.lats = [destination1[0], destination2[0]]
        #    self.lons = [destination1[1], destination2[1]]
        #else:
        #    self.lats = [destination2[0], destination1[0]]
        #    self.lons = [destination2[1], destination1[1]]           
        """
        probably also need to bound these between -180..180 and -90..90
        """
        """
        TRY ALTERNATIVE IMPLEMENTATION USING HAVERSINE FORUMLA
        http://www.movable-type.co.uk/scripts/latlong.html
        http://stackoverflow.com/questions/7222382/get-lat-long-given-current-point-distance-and-bearing
        """
        R = 6378137.0 #Radius of the Earth [m]
        brng = self.heading - 90
        brng = (brng + 360) % 360
        d = self.swath_width / 2.0 #Distance [m]
        lat1 = np.radians(self.lat)
        lon1 = np.radians(self.long)
        lat2 = np.arcsin( np.sin(lat1)*np.cos(d/R) +
             np.cos(lat1)*np.sin(d/R)*np.cos(brng))

        lon2 = lon1 + np.arctan2(np.sin(brng)*np.sin(d/R)*np.cos(lat1),
                             np.cos(d/R)-np.sin(lat1)*np.sin(lat2))
        # repeat with reverse Bearing
        brng = self.heading + 90
        brng = (brng + 360) % 360
        lat3 = np.arcsin( np.sin(lat1)*np.cos(d/R) +
             np.cos(lat1)*np.sin(d/R)*np.cos(brng))

        lon3 = lon1 + np.arctan2(np.sin(brng)*np.sin(d/R)*np.cos(lat1),
                             np.cos(d/R)-np.sin(lat1)*np.sin(lat2))
        # use as edges
        #self.lats = np.degrees([lat2, lat3])
        #self.lons =  np.degrees([lon2, lon3])


    def only_daytime_overpasses(self, sun):
        """
        Return only those which are visible -- this will take awhile i imagine
        """
        observer_on_the_ground = ephem.Observer()
        observer_on_the_ground.lat = np.radians(self.lat)
        observer_on_the_ground.long = np.radians(self.long)
        observer_on_the_ground.date = self.datetime

        """
        Improvment from
        http://stackoverflow.com/questions/15044521/
                javascript-or-python-how-do-i-figure-out-if-its-night-or-day

        think should work without timezones
        -- essentially next sunset should be before the next sunset

        For it to be the day time the sunset must be nearer than
        the next sunrise


        NOTE: towards the poles the sun never rises or sets
              depending on time of year. When this happens
              pyephem raises an error...
              These are some annoying edge cases to deal
              with

        """

        #import pdb; pdb.set_trace()
        try:
            next_sunrise= observer_on_the_ground.next_rising(sun).datetime()
            next_sunset = observer_on_the_ground.next_setting(sun).datetime()
        except (ephem.AlwaysUpError):
            # polar summer and sun never sets
            # so it must be daytime...
            self.daytime = True
            sun.compute(observer_on_the_ground)
            # can check with solarzenith angle
            solZenith = np.degrees(sun.alt)
            return None
        except (ephem.NeverUpError):
            # polar winter and sun never sets
            # so it must be nighttime
            self.daytime = False
            sun.compute(observer_on_the_ground)
            # can check with solarzenith angle
            solZenith = np.degrees(sun.alt)
            return None
        # if no exceptions it's not polar and sun does rise and set...
        if next_sunset < next_sunrise:
            self.daytime = True
        elif next_sunset > next_sunrise:
            self.daytime = False
        else:
            import pdb; pdb.set_trace() # something wrong
            self.daytime = -1 #error
        sun.compute(observer_on_the_ground)
        # can check with solarzenith angle
        solZenith = np.degrees(sun.alt)
        self.solarZenith = solZenith
        return None


def get_overpass(orbital, time):
    """
    Constructor function for making an overpass instance

    Function which returns a lat and lon for a time
    for the given pyephem orbital
    """
    orbital.compute(time)
    west = orbital.sublong < 0
    the_overpass = overpass( orbital.sublat / ephem.degree,
                             orbital.sublong / ephem.degree,
                             time, orbital.elevation, orbital.eclipsed,
                             west)
    return the_overpass

