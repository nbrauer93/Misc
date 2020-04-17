#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 22:50:20 2020

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
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



file = 'oisst-avhrr-v02r01.20180913.nc'
nc = Dataset(file, 'r')

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]-360

sst = nc.variables['sst'][0,0,:,:] #1,1,720, 1440 shape
anom = nc.variables['anom'][0,0,:,:]

lat2, lon2= np.meshgrid(lat,lon)




#%%
#Now let's plot


plt.figure(figsize=(10,10))
cmin = -4.; cmax = 4.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
xlim = np.array([-83,-60]); ylim = np.array([20,45])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates()
parallels = np.arange(20,46, step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-85, -59, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])

cs = m.contourf(lon2,lat2,anom.T,clevs,cmap='bwr',extend='both') 
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel(r'$^{o}$C',name='Calibri',size=20)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
cbar.set_ticks(clevs[::4])
cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)


  
#plt.title('Precipitation 8/27/2017',name='Calibri',weight='bold',size=20)
plt.title('Sea Surface Temperature Anomaly 13 Sept 2018',name='Calibri',size=20, y = 1.04)
plt.suptitle('1971-2000 climatology', size = 16, y= 0.91)
plt.show(block="False")