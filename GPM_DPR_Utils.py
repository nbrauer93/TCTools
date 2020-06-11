#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 13:31:26 2020

@author: noahbrauer
"""

import h5py 
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj
import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
  
import pyart

file = '2A.GPM.DPR.V8-20180723.20150620-S043651-E060924.007430.V06A.HDF5'
DPR = h5py.File(file, 'r')

lat = DPR['NS']['Latitude'][:,:]    
lon = DPR['NS']['Longitude'][:,:]
z = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:] #nscan x nray (7934,49)
dsd = DPR['NS']['SLV']['paramDSD'][:] #nscan x nray x nbin x DSDmoment (7934,49,176,2)
precip = DPR['NS']['SLV']['precipRateNearSurface'][:]
freezing = DPR['NS']['VER']['heightZeroDeg'][:]

zcolor = 'pyart_NWSRef'

def compute_xsect(raw_lat, raw_lon, start_lon, end_lon, ray_number):
    
    
    ind1 = np.where((lon[:,0]>=start_lon))
    ind2 = np.where((lon[:,0]<=end_lon))
    ind3 = np.intersect1d(ind1,ind2)
    
    x = 2.* 17 #48 degrees (from -17 to 17)
    re = 6378. #radius of the earth
    theta = -1*(x/2.) + (x/48.)*np.arange(0,49) #Split into equal degrees (from -17 to 17)
    theta2  = np.ones(theta.shape[0]+1)*np.nan #Define an empty array (NaNs) with same shape as ray dimension
    theta = theta - 0.70833333/2. #Shift degrees for plotting pruposes
    theta2[:-1] = theta #remove last array dimension in theta so python won't get angry about shape
    theta2[-1] = theta[-1] + 0.70833333
    theta = theta2*(np.pi/180.) #Convert from degrees to radians

    prh = np.ones((177,50))*np.nan #Define empty grid
    
    for i in range(prh.shape[0]): #Loop over each range gate
        for j in range(prh.shape[1]): #Loop over each scan
            a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #Orbit height of 407 km

            prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a)

    h2 = prh #where prh is the (range bins,ray)-space
    h3 =np.ones((177,50))*np.nan
    for i in range(h3.shape[1]):
        h3[:,i] = h2[::-1,i] #This reverses the vertical dimension so as indices increase, height increases
        
        
    ku = DPR['NS']['SLV']['zFactorCorrected'][ind3,:,:] #Read in ku-band reflectivity; nscan x nray (554,49,176)
    n0 = dsd[ind3,:,:,0] #Read in the number concentration
    d0 = dsd[ind3,:,:,1] #Read in the mean drop diameter  #Both have dimensions nscan x nray x nbin (554,49,176)
    zeroDeg = freezing[ind3,:]
    #Cut all parameters so they are at same ray as above
    
    ray = ray_number
    ku = ku[:,ray,:]
    n0 = n0[:,ray,:]
    d0 = d0[:,ray,:]
    zero_deg_isotherm = zeroDeg[:,ray]/1000 #Convert from meters to km
    
    lons = raw_lon[ind3,ray]
    lats = raw_lat[ind3,ray]
    
    #Choose a starting point, then calculate distance
    lat0 = lats[0]
    lon0 = lons[0]


    p = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0) #Define a projection and set starting lat an lon to same point as above  

    lat_3d = np.ones(ku.shape)*np.nan
    lon_3d = np.ones(ku.shape)*np.nan

    for i in range(ku.shape[0]):
        lat_3d[i,:] = lats[i]
        lon_3d[i,:] = lons[i]
        
    x,y = p(lon_3d,lat_3d) #Now convert degrees to distance (in km)
    R_gpm = np.sqrt(x**2 + y**2)*np.sign(x) #Keeps sign of number; converts to radial distance
    
    #Reverse range gate order for all parameters

    ku = ku[:,::-1]
    n0 = n0[:,::-1]
    d0 = d0[:,::-1]



    ku = np.ma.masked_where(ku<=12, ku) #Mask all the bad points in ku data
    y = np.ones([ku.shape[0], ku.shape[1]]) #Define an empty array
   
    
    #Define the appropriate range bins
    h4 = h3[:,ray] #This number MUST match the same ray being used
    for i in range(y.shape[1]):
        y[:,i] = h4[i]
        
    #Remove the values less than or equal to zero
    
    n0[n0<=0] = np.nan
    d0[d0<=0] = np.nan
    
    
    return R_gpm, y, d0, n0, ku, zero_deg_isotherm


data = compute_xsect(lat,lon, -90.5, -85.5, 34)
R_gpm = data[0]
y = data[1]
d0 = data[2]
n0 = data[3]
ku = data[4]
zero_deg_isotherm = data[5]



def plot_DPR (R_gpm, y, data_variables, zero_deg_isotherm, value_min, value_max, x_min, x_max, y_min, y_max, value_interval, label_font_size, title_font_size, colormap):
    
    
    plt.figure(figsize=(10,10))

    vmax = value_max
    vmin = value_min


    label_size = label_font_size
    title_size = title_font_size
    
    title_label = input("Enter the title label:")
    title = title_label
    
    xlabel = input("Enter the x-axis label:")
    x_label = xlabel
    
    ylabel = input("Enter the y-axis label")
    y_label = ylabel
    
    clabel = input("Enter the colorbar label:")
    cbar_label = clabel


    pm = plt.pcolormesh(R_gpm/1000., y, data_variables, cmap=colormap,vmin=vmin,vmax=vmax)
    pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k', label = r'$0^{o}$C isotherm')
    plt.xlabel(x_label, size = label_size)
    plt.ylabel(y_label, size = label_size)
    plt.title(title, size = title_size)
    plt.xlim(x_min,x_max)
    plt.ylim(y_min,y_max)
    plt.colorbar().set_label(label = cbar_label, size = label_size)
    plt.clim(vmin,vmax)
    plt.xticks(fontsize = label_size)
    plt.yticks(fontsize = label_size)

    return pm, pm2


plot = plot_DPR(R_gpm, y, ku, zero_deg_isotherm, 12, 60, 0, 250, 0, 15, 2.5, 24, 26, zcolor)

    
    



    
    

