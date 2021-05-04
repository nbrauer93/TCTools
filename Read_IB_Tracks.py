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
write_line = "storm_name,storm_date,storm_lat,storm_lon,sat_name,gpm_lat,gpm_lon,gpm_date,gpm_dist\n"
    
for i in range(len(time_numpy)):
    


    storm_name = name[i]
    storm_date = time_numpy[i]
    storm_lat = latitude[i]
    storm_lon= longitude[i]
    
    print(storm_name)
    print(storm_date)
    print(storm_lat)
    print(storm_lon)

    sat_name = "GPM"
    sat_swath = 245               # swath of DPR
    tle_list = "tle.GPM"          # name of TLE list file
    timewindow= 3               # Find all overpasses within +/- this many hours


    found_tle_file = "found_tle.txt"
    if ( os.path.exists(found_tle_file)): os.remove(found_tle_file)

    tlef= get_TLE(tle_list, storm_date, sat_name, found_tle_file)
    
    
    
    print(tlef)

    if ( os.path.exists(found_tle_file) ):
        print("Found TLE nearby date ", storm_date)

        rscat= Orbital(sat_name, tle_file=found_tle_file)

        passes= find_overpasstimes(rscat, storm_name, storm_lat, storm_lon, storm_date, timewindow, sat_name, sat_swath)
        
     
        
        

        for idx, lpass in enumerate(passes):
            print("Overpass number ",idx, " Now get the data according to these coordinates and times ", lpass)

            #print(passes)
            
        if len(passes) > 0:
            write_line += f"{storm_name},{storm_date},{storm_lat},{storm_lon},{passes[0][0]},{passes[0][1]},{passes[0][2]},{passes[0][3]},{passes[0][4]}\n"


o = open("output.txt","w")
o.write(write_line)
o.close()




#%%

#Now let's re-format the date-time indices from the output.txt file to the DPR filename format

#import wget


file_gpm = 'output.csv'


gpm_files = pd.read_csv(file_gpm)

datetime_storm = gpm_files['storm_date'][:]
datetime_gpm = gpm_files['gpm_date'][:]



#Loop through all GPM DPR dates/times and parse to match DPR file format
#### YYYYMMDD-SHHMMSS-EHHMMSS ##### 

dpr_file_times = np.asarray(datetime_gpm)

#Separate times and dates

date_only = [i[0:10] for i in dpr_file_times] 
    
times_only = [i[11:] for i in dpr_file_times]   


#Now remove hyphens and colons to conform to DPR filename format

#First hyphens from date_only:
date_gpm = []    

for i in date_only:

    date_gpm.append(i.replace('-',''))
    
#Ok now do the same for colons


time_gpm = []

for i in times_only:

    time_gpm.append(i.replace(':',''))    

 


#Also need to reformat time to get the orbit number to read into datetime

from datetime import datetime



utc_time = [datetime.strptime(date_gpm[i]+time_gpm[i][0:4],"%Y%m%d%H%M") for i in range(len(date_gpm))]

#Now lets extract the orbit number


orbits = []

write_orbit_no = "orbit_no"


for i in range(len(utc_time)):
    
    orbit_no = rscat.get_orbit_number(utc_time[i])
    print(orbit_no)
    
    orbits.append(orbit_no)
    
    #Write to a text file
    write_orbit_no += f"{orbit_no}\n"

    
w = open("orbit_no.txt", 'w')
w.write(write_orbit_no)
w.close()














