#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 20:40:12 2022

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset, num2date, MFDataset
import netCDF4 as nc
import glob
import imageio
import os


files = glob.glob('*.nc')

#Let's parse the end of the file name, then sort to grab indices

times_parsed = []
times_string = []

for i in range(len(files)):  
    
    first_split = files[i].split("_")[3]
    #(files[i])
    #print(first_split)
    second_split = first_split.split(".")[0]
    #print(second_split)
    
    
    #Remove first zero from string:
        
    third_split = second_split[0:4]
    print(third_split)
    times_parsed.append(third_split)
    
    times_used = third_split[0:4]
    #print(times_used)
    times_string.append(times_used)
   

#Now sort times

times_ordered_index = np.argsort(times_parsed)    


times_ordered = [files[i] for i in times_ordered_index]




#%%


#order time

time_array = np.array(times_string)
time_array_sort = [time_array[i] for i in times_ordered_index]


    
#%%
#Open one file to retrieve lat-lon

file = 'KTBW_N0Q_20220928_150000.nc'

nc = Dataset(file, 'r')

latitude = nc.variables['lat'][:]
longitude = nc.variables['lon'][:]

lat2,lon2 = np.meshgrid(latitude, longitude)




for file in range(len(times_ordered)):
    
    data = Dataset(times_ordered[file], 'r')
    
    zh = data.variables['bref'][:]
    
    zh[zh<0] = np.nan
    
    #Now plot
    plt.figure(figsize=(20,20))
    
    cmin = 0; cmax = 60; cint = 2; clevs = np.round(np.arange(cmin,cmax,cint),2)
    xlim = np.array([-85,-78.5]); ylim = np.array([24,30])
   
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
    m.drawstates(); m.drawcountries(); m.drawcoastlines()
    m.drawcounties()
    cs = m.contourf(lon2,lat2,zh[0,:,:].T, clevs, cmap = 'turbo', extend = 'both')
    
    
    x_loc = -80.75
    y_loc = 23.5
    label = '@NOAABrauer'
    plt.text(x_loc, y_loc, label, size = 24)
    
   

    cbar = plt.colorbar(fraction=0.046)
    cbar.ax.tick_params(labelsize = 26)
    cbar.set_label(label = '[dBZ]',size = 26)
        
    plt.title(r'KTBW 09/28 ' + time_array_sort[file] + ' UTC' + ' $0.5^{o} Z_{H}$', size = 40)
    plt.savefig(time_array_sort[file]+'.png')
    plt.show()
    