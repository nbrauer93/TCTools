#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:43:13 2020

@author: noahbrauer
"""
import numpy as np
import matplotlib.pyplot as plt

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


def draw_vert_line(lat_min,lat_max,longitude):
    
    latitude_vert = np.arange(lat_min,lat_max,step = 1)
    longitude_vert = np.full((lat_max-lat_min),longitude)
    
    return longitude_vert,latitude_vert

def draw_horiz_line(lon_min,lon_max,latitude):
    
    longitude_horiz = np.arange(lon_min,lon_max,step = 1)
    latitude_horiz = np.full((lon_max-lon_min),latitude)
    
    return longitude_horiz,latitude_horiz


def positive_slope(lat_min,lat_max,lon_min,lon_max):
    
    lat = [lat_min,lat_max]
    lon = [lon_min,lon_max]
    
    return lon,lat


def negative_slope(lat_min,lat_max,lon_min,lon_max):
    
    lat = [lat_min,lat_max]
    lon = [lon_max,lon_min]
    
    return lon,lat


#%%

plt.figure(figsize=(14,14))

xlim = np.array([-180,180]); ylim = np.array([-40,40])
   
m = Basemap(projection='merc',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(); m.drawcountries()

parallels = np.arange(-50,50, step = 10)
m.drawparallels(parallels, labels = [True, True, False, False])

meridians = np.arange(-180, 180, step = 20)
m.drawmeridians(meridians, labels = [False, False, False, True])

linewidth = 4  
font_size = 17
name = 'Calibri'
    
vertical_pac = draw_vert_line(-40,40,100)
x,y = m(vertical_pac[0],vertical_pac[1])
m.plot(x,y, color = 'k', linewidth = linewidth)   


horiz_pac = draw_horiz_line(105,180,0)
x,y = m(horiz_pac[0],horiz_pac[1])
m.plot(x,y, color = 'k', linewidth = linewidth)

horiz_indian = draw_horiz_line(42,105,0)
x,y = m(horiz_indian[0],horiz_indian[1])
m.plot(x,y,color = 'k', linewidth = linewidth)

horiz_pac = draw_horiz_line(-180,-80,0)
x,y = m(horiz_pac[0],horiz_pac[1])
m.plot(x,y, color = 'k', linewidth = linewidth)

vert_pac = draw_vert_line(0,40,-179)
x,y = m(vert_pac[0],vert_pac[1])
m.plot(x,y,color = 'k', linewidth = linewidth)


lon = negative_slope(10,40,-120,-83)[0]
lat = negative_slope(10,40,-120,-83)[1]
x,y = m(lon,lat)
m.plot(x,y, color = 'k', linewidth = linewidth)


lon = positive_slope(-33,0,28,42)[0]
lat = positive_slope(-33,0,28,42)[1]
x,y = m(lon,lat)
m.plot(x,y,color = 'k', linewidth = linewidth)

vertical_indian = draw_vert_line(0,30,42)
x,y = m(vertical_indian[0],vertical_indian[1])
m.plot(x,y,color = 'k', linewidth = linewidth)


###Add Text
x2star,y2star = m(-155,32)
plt.text(x2star,y2star,'ECPAC' , color = 'k', fontsize = font_size)

x3star,y3star = m(-150,-24)
plt.text(x3star,y3star,'SPAC' , color = 'k', fontsize = font_size)

x4star,y4star = m(-54,26)
plt.text(x4star,y4star,'ATL' , color = 'k', fontsize = font_size)

x5star,y5star = m(57.5,5)
plt.text(x5star,y5star,'NIND' , color = 'k', fontsize = font_size)

x6star,y6star = m(69,-23)
plt.text(x6star,y6star,'SIND' , color = 'k', fontsize = font_size)

x7star,y7star = m(160,-31)
plt.text(x7star,y7star,'SPAC' , color = 'k', fontsize = font_size)

x8star,y8star = m(145,22)
plt.text(x8star,y8star,'NWPAC' , color = 'k', fontsize = font_size)









