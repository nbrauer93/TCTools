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

import pandas as pd


#Now open the CSV file



ibtracks_data = pd.read_csv('ibtracs.ALL.list.v04r00.csv')

#Extract time attributes


season = ibtracks_data['SEASON'][1:]
season = season.copy()

season_array = pd.Series.to_numpy(season, dtype = 'float')




season_index = np.where(season_array>2013)[0]


season_GPM = season_array[season_index]


#Extract lats and lons:

lat = ibtracks_data['LAT'][season_index]
lon = ibtracks_data['LON'][season_index]


latitude = pd.Series.to_numpy(lat,dtype = 'float')
longitude = pd.Series.to_numpy(lon,dtype = 'float')


#Extract basin and name

basin = ibtracks_data['BASIN'][season_index]
name = ibtracks_data['NAME'][season_index]


name = name.copy()

#NaN out all the unnamed storms

name[name == "NOT_NAMED"] = np.nan

name = pd.Series.to_numpy(name, dtype = 'str')

#Extract exact times


time = ibtracks_data['ISO_TIME'][season_index]
time_numpy = pd.Series.to_numpy(time, dtype = 'str')


#%%
###Now lets loop through each time and find all the GPM files that match up:
    
    
for i in range(len(time_numpy)):
    
    
    
    if __name__ == '__main__':

        storm_name = name[i]
        storm_date = time_numpy[i]
        storm_lat = latitude[i]
        storm_lon= longitude[i]
        
        print(storm_name)
        print(storm_date)
        print(storm_lat)
        print(storm_lon)

        sat_name = "GPM"
        sat_swath = 800               # swath of GMI
        tle_list = "tle.GPM"          # name of TLE list file
        timewindow= 24                # Find all overpasses within +/- this many hours


        found_tle_file = "found_tle.txt"
        if ( os.path.exists(found_tle_file) ): os.remove(found_tle_file)

        tlef= get_TLE(tle_list, storm_date, sat_name, found_tle_file)
        
        
        print(tlef)

        if ( os.path.exists(found_tle_file) ):
            print("Found TLE nearby date ", storm_date)

            rscat= Orbital(sat_name, tle_file=found_tle_file)

            passes= find_overpasstimes(rscat, storm_name[i], storm_lat[i], storm_lon[i], storm_date[i], timewindow, sat_name, sat_swath)

            for idx, lpass in enumerate(passes):
                print("Overpass number ",idx, " Now get the data according to these coordinates and times ", lpass)

                #print(passes)

    if ( os.path.exists(found_tle_file) ): os.remove(found_tle_file)
    
    gpm_tle = []
    gpm_tle.append(found_tle_file)

        
