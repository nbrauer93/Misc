#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 12:08:46 2023

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
import pyart
import cartopy.crs as ccrs
import glob
#from tqdm import tqdm
from scipy import stats

radar_file = 'KMPX20230223_140917_V06'


#import numpy as np

#from pyart.io import read
#from pyart.core import antenna_to_cartesian



def quasi_vertical_profile(filename, fields=None, gatefilter=None):

     radar = pyart.io.read(filename)

     if fields is None:
         fields = radar.fields

    #if gatefilter is not None:


     desired_angle = 20.0
     index = abs(radar.fixed_angle['data'] - desired_angle).argmin()
     #print(radar.fixed_angle['data'])
     #print(radar.elevation['data'][-1])

     qvp = {}

     for field in fields:
         this_field = radar.get_field(index, field).mean(axis = 0)
         qvp.update({field:this_field})

     qvp.update({'range': radar.range['data'], 'time': radar.time})
     x,y,z = pyart.core.antenna_to_cartesian(qvp['range']/1000.0, 0.0,
                                 radar.fixed_angle['data'][index])
     qvp.update({'height': z})
     #del radar
     return qvp

qvp = quasi_vertical_profile(radar_file)
#print(qvp['reflectivity'])

plt.plot(qvp['differential_reflectivity'],qvp['height']/1000.0)
plt.xlabel('Mean Differential Reflectivity (dB)')
plt.ylabel('Altitude')
plt.title('Quasi-Vertical Profile 2154 UTC')
#plt.show()


#Now do this for multiple files



files = glob.glob('KMPX*')

#Loop through file list to get times; Parse out of file radar_name
#Sort file names by time first:

order_index = np.argsort(files)
files_ordered = [files[i] for i in order_index]

times = []
qvps_z = np.ones((136,1832))*np.nan
qvps_zdr = np.ones((136,1832))*np.nan
qvps_rhohv = np.ones((136,1832))*np.nan

for i in range(len(files_ordered)):

    qvps = quasi_vertical_profile(files_ordered[i])

    qvps_z[i,:] = qvps['reflectivity']
    qvps_zdr[i,:] = qvps['differential_reflectivity']
    qvps_rhohv[i,:] = qvps['cross_correlation_ratio']

    parsed_time = files_ordered[i].split("_")[1]
    print(parsed_time)
    times.append(parsed_time)

#print(qvps_z.shape) # 107 x 1832
#print(np.nanmax(qvps_z))  # 44.495 dBZ

times = np.asarray(times)



#%%


#Now plot:

font_size = 16    
title_size = 20

plt.figure(figsize=(10,10))
    
qvps_z[qvps_rhohv<0.8] = np.nan
qvps_zdr[qvps_rhohv<0.8] = np.nan
qvps_rhohv[qvps_rhohv<0.8] = np.nan

plt.pcolormesh(times,qvp['height']/1000.0,qvps_z.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size= font_size)
plt.ylabel('Height (km)', size = font_size)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = '[dBZ]',size = font_size)
#plt.colorbar()
plt.clim(-15,30)
plt.ylim(0,10)

plt.title(r'KMPX Quasi-Vertical Profiles 02/22-23 $Z_{H}$', size = title_size)
x = [0,17,34,51,68,85,102,119,136]
labels = np.array(['2000','2230','0100','0330','0600','0830','1100','1330','1600'])
plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()


plt.figure(figsize=(10,10))

plt.pcolormesh(times,qvp['height']/1000.0,qvps_zdr.T, cmap = 'pyart_RefDiff')
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('Height (km)', size = font_size)
plt.ylim(0,10)
plt.clim(-0.5,3)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = '[dB]',size = font_size)
plt.title(r'KMPX Quasi-Vertical Profiles 02/22-23 $Z_{DR}$', size = title_size)

plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()

plt.figure(figsize=(10,10))

plt.pcolormesh(times,qvp['height']/1000.0,qvps_rhohv.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('Height (km)', size = font_size)
plt.ylim(0,10)
plt.clim(0.9,1.)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = r'[$\rho_{HV}$]',size = font_size)
plt.title(r'MPX Quasi-Vertical Profiles 02/22-23 $\rho_{HV}$', size = title_size)

plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()