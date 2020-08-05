#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:29:22 2020

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import netCDF4 as nc
import numpy as np
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator




import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
import pyart
### File includes September 1-30 0000 UTC (1989-2019)  850 mb and 200 mb winds 


file = 'sept_global_shear.nc'

#Read in the file

nc = Dataset(file, 'r')

#Read in file attributes

longitude = nc.variables['longitude'][:]
latitude = nc.variables['latitude'][:]


level = nc.variables['level'][:]

u_200 = nc.variables['u'][:,0,:,:]*1.94384
v_200 = nc.variables['v'][:,0,:,:]*1.94384

u_850 = nc.variables['u'][:,1,:,:]*1.94384
v_850 = nc.variables['v'][:,1,:,:]*1.94384

#Create grid for plotting

lat2,lon2 = np.meshgrid(latitude,longitude)


era = {}

time = nc.variables['time'][:]
timeUnits = nc.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
era['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
era['day'] = np.asarray([d.day for d in era['date']])
era['month'] = np.asarray([d.month for d in era['date']])
era['year'] = np.asarray([d.year for d in era['date']])


#Calculate the magnitude of the wind

def wind_magnitude(u,v):
    
    magnitude = np.sqrt((u**2)+(v**2))
    
    return magnitude

wind_850 = wind_magnitude(u_850,v_850)
wind_200 = wind_magnitude(u_200,v_200)


#Compute wind shear: The vector difference between the 200 mb wind and the 850 mb wind

def wind_shear(wind_upper,wind_lower):
    shear = wind_upper-wind_lower
    
    return shear


shear_mag = wind_shear(wind_200,wind_850)



#Now compute the mean and standard deviations of wind shear magnitude from 1989-2019 for September

shear_mean_period = np.nanmean(shear_mag, axis = 0)
shear_std_period = np.nanstd(shear_mag, axis = 0)

#Index only 2014-2019 years (since the launch of GPM)

gpm_index = np.where(era['year']>2013)[0]

shear_gpm = shear_mag[gpm_index,:,:]

#Daily standardized shear anomalies: 0 index is September 1, 2014
time_index = 0

z_score_sept1 = (shear_gpm[time_index,:,:] - shear_mean_period)/shear_std_period

#%%


title_name = '850-200 mb Shear Anomaly '
time = '9/1/2014 0000 UTC'  
title_font_size = 22 


cmin = -4.5; cmax = 4.5; cint = 0.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
   
plt.figure(figsize=(20,10))


xlim = np.array([-180,180]); ylim = np.array([-50,50])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(); m.drawcountries()

parallels = np.arange(-50,50, step = 10)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-180, 180, step = 20)
m.drawmeridians(meridians[::2], labels = [False, False, False, True])


cs = m.contourf(lon2,lat2,z_score_sept1.T, clevs, cmap = 'bwr', extend = 'both')
    

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel(r'[$\sigma$]',size=title_font_size)
plt.title(str(title_name) + str(time),name='Calibri',size=title_font_size)
plot = plt.show(block=False) 




shear_mean_period[shear_mean_period<5] = np.nan

title_name = 'September Mean 850-200 mb Shear (1989-2019) '
title_font_size = 22 


cmin = 5; cmax = 70; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
   
plt.figure(figsize=(20,10))


xlim = np.array([-180,180]); ylim = np.array([-50,50])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(); m.drawcountries()

parallels = np.arange(-50,50, step = 20)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-180, 180, step = 20)
m.drawmeridians(meridians[::2], labels = [False, False, False, True])


cs = m.contourf(lon2,lat2,shear_mean_period.T, clevs, cmap = 'YlOrRd', extend = 'both')
    

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[knots]',size=title_font_size)
plt.title(str(title_name),name='Calibri',size=title_font_size)
plot = plt.show(block=False) 





