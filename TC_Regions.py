#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:33:01 2021

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Import file

file = 'output.csv'

data = pd.read_csv(file)

#Read in lat and lon

gpm_lat = data['gpm_lat'][:]
gpm_lon = data['gpm_lon'][:]

gpm_date = data['gpm_date'][:]

#Change to numpy arrays

lon = pd.Series.to_numpy(gpm_lon)
lat = pd.Series.to_numpy(gpm_lat)



#Determine booleans for each region based off lat-lon criteria

epac = np.logical_and(lat>0,np.logical_and(lon>-180,lon<-100))
spac = np.logical_and(lat<0, np.logical_or(lon<-60, lon>100))
atl = np.logical_and(lon>-100,lon<20)
nind = np.logical_and(lat>0,np.logical_and(lon>40,lon<100))
sind = np.logical_and(lat<0,np.logical_and(lon>25,lon<100))
nwpac = np.logical_and(lat>0,np.logical_and(lon>100,lon<180))


#Change booleans to integers

epac = epac.copy().astype(int)
spac = spac.copy().astype(int)
atl = atl.copy().astype(int)
nind = nind.copy().astype(int)
sind = sind.copy().astype(int)
nwpac = sind.copy().astype(int)


'''

EPAC = 1
SPAC = 2
ATL = 3
NIND = 4
SIND = 5
NWPAC = 6

'''

#Now loop through each list; If value = 1, assign to value corresponding to each region


regions = []


for i in range(len(epac)):
    
    if epac[i] == 1:
        
        value = 1
        regions.append(value)
        
    elif spac[i] == 1:
        
        value = 2
        regions.append(value)
        
        
    elif atl[i] == 1:
        
        value = 3
        regions.append(value)
        
        
    elif nind[i] == 1:
        
        value = 4
        regions.append(value)
        
    elif sind[i] == 1:
        
        value= 5
        regions.append(value)
        
    elif nwpac[i] == 1:
        
        value = 6
        regions.append(value)
        
    else: 
        
        value = np.nan
        regions.append(value)




        
        
  
    
    
