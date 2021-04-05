#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 15:23:27 2021

@author: noahbrauer
"""


#Import all libraries

import matplotlib.pyplot as plt
import numpy as np
import datetime
import re
import sys, os
from pyorbital.orbital import Orbital
from math import sin, cos, sqrt, atan2, radians

from overpasstimeswithlist import distance_km,get_TLE,find_overpasstimes

from mpl_toolkits.basemap import Basemap

import csv


#Now open the CSV file


with open('ibtracs.ALL.list.v04r00.csv', newline = '') as csvfile:
    
    ibtracs = csv.reader(csvfile, delimiter = ' ', quotechar = '|')
    
    for row in ibtracs:
        
        print(','.join(row))






#%%


if __name__ == '__main__':

  storm_name = "SALLY"
  storm_date = '2020-09-15 00:00:00'
  storm_lat = 30.0
  storm_lon= -99.0

  sat_name = "GPM"
  sat_swath = 800               # swath of GMI
  tle_list = "tle.GPM"          # name of TLE list file
  timewindow= 24                # Find all overpasses within +/- this many hours


  found_tle_file = "found_tle.txt"
  if ( os.path.exists(found_tle_file) ): os.remove(found_tle_file)

  tlef= get_TLE(tle_list, storm_date, sat_name, found_tle_file)

  if ( os.path.exists(found_tle_file) ):
    print("Found TLE nearby date ", storm_date)

    rscat= Orbital(sat_name, tle_file=found_tle_file)

    passes= find_overpasstimes(rscat, storm_name, storm_lat, storm_lon, storm_date, timewindow, sat_name, sat_swath)

    for idx, lpass in enumerate(passes):
      print("Overpass number ",idx, " Now get the data according to these coordinates and times ", lpass)

    #print(passes)

  if ( os.path.exists(found_tle_file) ): os.remove(found_tle_file)