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
    r"""
    This function calculates the magnitude of the  wind from both u and v components. Outputs values in m/s.
    
    
    Parameters:
    -----------
    u: float
        U-wind in m/s
    v: float
        V-wind in m/s
    
    """
    
    magnitude = np.sqrt((u**2)+(v**2))
    
    return magnitude


#Compute magnitude of shear vector

def wind_shear(u_upper,u_lower, v_upper, v_lower):
    r"""
    This function computes the magnitude of the vertical wind shear, the magnitude of the lower wind, and the magnitude of the upper wind.
    
    Parameters:
    -----------
    u_upper: float
        U-wind at the upper level in m/s
    u_lower: float
        U-wind at the lower level in m/s
    v_upper: float
        V-wind at the upper level in m/s
    v_lower: float
        V-wind at the lower level in m/s
    
    
    """
    
    magnitude_upper =  np.sqrt((u_upper**2)+(v_upper**2))
    magnitude_lower = np.sqrt((u_lower**2)+(v_lower**2))
    
    shear = magnitude_upper - magnitude_lower

    
    return shear, magnitude_upper, magnitude_lower


#Compute  u and v components of wind shear vector

def wind_shear_comp_vectors (u_upper, u_lower, v_upper, v_lower):
    r"""
    
    """
    
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
    r"""
    This function transforms lat-lon distance into metric distance (in kilometers)
    
    Parameters:
    -----------
    lat_input2: float
        Latitude coordinate of point two
    lat_input1: float
        Latitude coordinate of point one
    lon_input2: float
        Longitude coordinate of point two
    lon_input1: float
        Longitude coordinate of point one
    
    
    """
    
    distlat = (lat_input2*np.cos(lat_input2)-lat_input1*np.cos(lat_input1))*111
    distlon = (lon_input2 - lon_input1)*111
    totaldist = np.sqrt((distlat)**2 + (distlon)**2)
    
    return(totaldist)
    
def perp_vector(u_shear, v_shear):
    
    r"""
    This function calculates the vector orthogonal to the shear vector, the slope of the perpendictular vector, and each component of the vector.
    
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
    

def shear_rel_quads(shear_vector, perp_vector, radar_data):
    
    r"""
    This function computes shear-relative quadrants of the tropical cyclone. The coordinate for the HURDAT2 Best Track location is the origin (0,0). 
    
    Parameters:
    -----------
    shear_vector: float (2D array)
        Shear vector
    coord: float (2D array)
        Coordinate of the center of the TC from HURDAT2
    perp_vector: float (2D array)
        Vector that is perpendicular to the shear vector
    radar_data: float (3D array; latxlonxtime at a constant elevation/altitude level)
    
    """
    
    #Define the origin as the intersection of the shear vector and the perpendicular vector; also the coordinate of the HURDAT2
    
    coord = np.array([0,0]) 
    
    radar_data = np.ones((radar_data.shape[0], radar_data.shape[1], radar_data.shape[2]))*np.nan
    downshear_left = np.ones((radar_data.shape[0], radar_data.shape[1], radar_data.shape[2]))*np.nan
    downshear_right = np.ones((radar_data.shape[0], radar_data.shape[1], radar_data.shape[2]))*np.nan
    upshear_left = np.ones((radar_data.shape[0], radar_data.shape[1], radar_data.shape[2]))*np.nan
    upshear_right = np.ones((radar_data.shape[0], radar_data.shape[1], radar_data.shape[2]))*np.nan
    
    #downshear_right[(perp_vector[0] > coord) & (shear_vector[1] > coord)] = radar_data[(perp_vector[0] > coord) & (shear_vector[1] > coord)]
    
    for i in range(radar_data.shape[0]):
        for j in range(radar_data.shape[1]):
            for k in range(radar_data.shape[2]):
                
 
                if perp_vector[0] and shear_vector[1] > coord:
                    downshear_right[i,j,k] = radar_data[i,j,k]
    
                elif perp_vector[0] and shear_vector[1]< coord:
                    upshear_left[i,j,k] = radar_data[i,j,k]
        
                elif perp_vector[0]>coord and shear_vector[1]<coord:
                    upshear_right[i,j,k] = radar_data[i,j,k]
    
                elif perp_vector[0]<coord and shear_vector[1]>coord:
                    downshear_left[i,j,k] = radar_data[i,j,k]
        
    return downshear_right, downshear_left, upshear_right, upshear_left




    
def annuli(coord, eyewall_radius, inner_band_radius, outer_band_radius):
    r"""
    Determines annuli (eyewall, inner, outer); From Kumjian et. al (2017); 15 km is eyewall, 60 km inner rainbands, 110 km is outer rainbands
    
    Parameters:
    -----------
    coord: float (2D lat-lon array)
        Coordinate of the center of the TC from HURDAT2
    eyewall_radius: float
        The user-set radius for the eyewall annulus
    inner_band_radius: float
        
    """
    
    
    #Convert km to lat-lon grid
    
    eyewall_regrid = eyewall_radius/111
    inner_regrid = inner_band_radius/111
    outer_regrid = outer_band_radius/111

    
    #Define the annuli for each

    
    

    