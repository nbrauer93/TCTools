#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 20:08:00 2020

@author: noahbrauer
"""

import numpy as np
import matplotlib.pyplot as plt



#Compute wind magnitude

def wind_magnitude(u,v):
    magnitude = np.sqrt((u**2)+(v**2))
    
    return magnitude


#Compute magnitude of shear vector

def wind_shear(u_upper,u_lower, v_upper, v_lower):
    
    magnitude_upper =  np.sqrt((u_upper**2)+(v_upper**2))
    magnitude_lower = np.sqrt((u_lower**2)+(v_lower**2))
    
    shear = magnitude_upper - magnitude_lower

    
    return shear, magnitude_upper, magnitude_lower


#Compute  u and v components of wind shear vector

def wind_shear_comp_vectors (u_upper, u_lower, v_upper, v_lower):
    
    u_shear = u_upper - u_lower
    v_shear = v_upper - v_lower
    shear_vector = np.array([u_shear,v_shear])
    
    return u_shear, v_shear, shear_vector



#Compute potential temperature (in degrees Kelvin)

def potential_temperature(temp,pressure):
    r"""
    This function calculates potential temperature and outputs value in Kelvin.
    
    Parameters:
    -----------
    temp : float
        Temp in Kelvin
    pressure : float
        Pressure in hPa
    """
    theta = temp*(1000/pressure)**0.286
    
    return theta

#Convert latitude/longitude to distance 

def latlon_to_dist(lat_input2,lat_input1,lon_input2,lon_input1):
    distlat = (lat_input2*np.cos(lat_input2)-lat_input1*np.cos(lat_input1))*111
    distlon = (lon_input2 - lon_input1)*111
    totaldist = np.sqrt((distlat)**2 + (distlon)**2)
    
    return(totaldist)
    
def perp_vector(u_shear, v_shear):
    
    r"""
    This function calculates the vector orthogonal to the shear vector
    
    Parameters:
    -----------
    u_shear: float
        U-component of the shear vector in m/s
    v_shear: float
        V-component of the shear vector in m/s
    
    """
    
    
    #Compute negative reciprocal of shear vector slope; this is the slope of the perpendicular line:
    
    slope = v_shear/u_shear
    perp_line_slope = -u_shear/v_shear
    perp_line_i = v_shear
    perp_line_j = -u_shear    
    perp_vector = np.array([v_shear, -u_shear])
    return perp_line_slope, perp_line_i, perp_line_j, perp_vector
   
    
    
#Determine shear-relative quadrants based off HURDAT2 data; coord should be a lat-lon coordinate, shear is a 2D vector 
    

def shear_rel_quads(shear_vector, coord, perp_vector):
    
    r"""
    This function computes shear-relative quadrants of the tropical cyclone
    
    Parameters:
    -----------
    shear_vector: float (2D array)
        Shear vector
    coord: float (2D array)
        Coordinate of the center of the TC from HURDAT2
    perp_vector: float (2D array)
        Vector that is perpendicular to the shear vector
    
    """
    
    #Define the origin as the intersection of the shear vector and the perpendicular vector
    
    intersection = np.intersect2d(shear_vector,perp_vector)
    
    
    
    
    
    




#Determine annuli (eyewall, inner, outer); From Kumjian et. al (2017); 15 km is eyewall, 60 km inner rainbands, 110 km is outer rainbands
    
def annuli(coord, eyewall_radius, inner_band_radius, outer_band_radius):
    
    
    #Convert km to lat-lon grid
    
    eyewall_regrid = eyewall_radius/111
    inner_regrid = inner_band_radius/111
    outer_regrid = outer_band_radius/111

    
    #Define the annuli for each

    
    

    