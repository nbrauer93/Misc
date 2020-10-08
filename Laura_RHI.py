#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 13:53:14 2020

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

from netCDF4 import Dataset, num2date, MFDataset





import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

import pyart

#Read in file

file = 'ncswp_SR1_20200827_011704.017_v125_s000_160.0_RHI_.nc'

#Import file

nc = Dataset(file, 'r')


#Extract attributes in file

zh = nc.variables['DBZ'][:]
zdr = nc.variables['ZDR'][:]
kdp = nc.variables['KDP'][:]
rho_hv = nc.variables['RHOHV'][:]

zdr_calibration = nc.variables['zdr_correction_db'][:]


elevation = nc.variables['Elevation'][:]

#Now calculate height (height) and distance (x) from the radar using Doviak and Zrnic:


deg_to_rad = np.pi/180. #Convert from degrees to radians

range_gate = np.zeros_like(zh[0,:])
range_gate = np.arange(0,len(range_gate))
range_gate = nc.variables['Range_to_First_Cell'].getValue() + range_gate*nc.variables['Cell_Spacing'].getValue()
a_e = (8494.66667)*1000


height = np.ones(zh.shape)*np.nan
s = np.ones(zh.shape)*np.nan  

for i in range(0,len(elevation)):

    a = range_gate**2 + a_e**2 + (2.0*range_gate*a_e*np.sin(elevation[i]*deg_to_rad))
    height[i,:] =((a**0.5)-a_e)+2
    s[i,:] = a_e*np.arcsin((range_gate*np.cos(elevation[i]*deg_to_rad))/(a_e+height[i,:]))




#Average over each column
    
    
s_avg = np.nanmean(s, axis = 0)    


#Now let's plot        
    #%%
label_size = 20    
tick_label_declutter = 200

height_1d = height[:,1]
#Round s and h values to the nearest decimal place:

#s_rounded = np.around(s_avg,decimals = 1)

    
fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height, zh, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(5,70,step = 2.5))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylim(0,14000)

#Need to plot x-ticks every 25 km from the radar
#Loop through s_avg to determine which value in array is divislbe by 25 (as an integer)
#If yes, store to array called xticks_25km


xticks_25km = []


#for i in range(len(s_avg)):
    
    
    
    
    
    
    
    





plt.xticks(s_avg[::tick_label_declutter]/1000,s_avg[::tick_label_declutter]/1000, size = label_size)
ax.set_xticklabels(s_avg[::tick_label_declutter]/1000)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

plt.yticks(height_1d[::tick_label_declutter], height_1d[::tick_label_declutter])
ax.set_yticklabels(height_1d[::tick_label_declutter])
plt.title(r'SMART-R2 $160^{o}$ Azimuth $Z_{H}$ 8/27 0117 UTC', size = label_size)
plt.show()




fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height, zdr-2.5, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(0,6,step = 0.25))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dB]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylim(0,14000)
plt.xticks(s_avg[::tick_label_declutter]/1000,s_avg[::tick_label_declutter]/1000, size = label_size)
ax.set_xticklabels(s_avg[::tick_label_declutter]/1000)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))


plt.yticks(height_1d[::tick_label_declutter], height_1d[::tick_label_declutter])
ax.set_yticklabels(height_1d[::tick_label_declutter])
plt.title(r'SMART-R2 $160^{o}$ Azimuth $Z_{DR}$ 8/27 0117 UTC', size = label_size)
plt.show()





    
fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height, rho_hv, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(0.9,1.05,step = 0.001))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[$\rho_{hv}$]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylim(0,14000)
plt.xticks(s_avg[::tick_label_declutter]/1000,s_avg[::tick_label_declutter]/1000, size = label_size)
ax.set_xticklabels(s_avg[::tick_label_declutter]/1000)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))


plt.yticks(height_1d[::tick_label_declutter], height_1d[::tick_label_declutter])
ax.set_yticklabels(height_1d[::tick_label_declutter])
plt.show()





print(np.nanmax(zh))
print(np.nanmax(zdr))
print(np.nanmax(rho_hv))
print(np.nanmax(kdp))












 



