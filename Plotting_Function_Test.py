#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:16:37 2020

@author: noahbrauer
"""

from GPM_DPR_Utils import compute_xsect, plot_DPR

import h5py 
import numpy as np
import matplotlib.pyplot as plt
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

data = compute_xsect(lat,lon, -90.5, -85.5, 34)
R_gpm = data[0]
y = data[1]
d0 = data[2]
n0 = data[3]
ku = data[4]
zero_deg_isotherm = data[5]

plot = plot_DPR(R_gpm, y, ku, zero_deg_isotherm, 12, 60, 0, 250, 0, 15, 2.5, 24, 26, zcolor)