#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 23:45:53 2022

@author: noahbrauer
"""


import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import pyart



file = 'RAXPOL-20220929-045523-A64.0-Z.nc'

#Read in file

nc = Dataset(file, 'r')

#Extract attributes:
    
z = nc.variables['Intensity'][:]
azimuth = nc.variables['Azimuth'][:]
elevation = nc.variables['Elevation'][:]

    
#Now calculate beam height and distance from the radar using Doviak and Zrnic



deg_to_rad = np.pi/180. #Convert from degrees to radians

range_to_first_cell = 0.
gate_width = nc.variables['GateWidth'][:]

range_gate = np.zeros_like(z[0,:])
range_gate = np.arange(0,len(range_gate))




range_gate_distance = range_gate.copy()
range_gate_distance = []

for i in range(len(range_gate)):
    
    dist = range_gate[i]*30
    range_gate_distance.append(dist)
    
    
range_gate_distance = np.asarray(range_gate_distance)    



#Effective radius 
a_e = (8494.66667)*1000 
    
#Now compute height/distance from the radar

height = np.ones(z.shape)*np.nan
s = np.ones(z.shape)*np.nan  

for i in range(0,len(elevation)):

    a = range_gate**2 + a_e**2 + (2.0*range_gate_distance*a_e*np.sin(elevation[i]*deg_to_rad))
    height[i,:] =((a**0.5)-a_e)+2
    s[i,:] = a_e*np.arcsin((range_gate_distance*np.cos(elevation[i]*deg_to_rad))/(a_e+height[i,:]))
    
#NaN out bad z values


z[z<0.1] = np.nan
    
#Now let's plot les data   
#%%    

label_size = 24
tick_label_declutter = 200



    
fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, z, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(5,70,step = 2.5))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)

xlabels = np.arange(-50,50,10)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)
plt.xlim(-50,50)
plt.xticks(xlabels, size = label_size)
plt.yticks(ylabels, size = label_size)

plt.title(r'RaXPol $64^{o}$ Azimuth $Z_{H}$ 9/29 0455 UTC', size = label_size)
plt.show()    
    
    
    
    
    
    
    
    