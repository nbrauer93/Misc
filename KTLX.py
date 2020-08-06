#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 17:08:10 2020

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

import pyart
from mpl_toolkits.basemap import Basemap

file = 'KTLX_N0R_20200712_235100.nc'

nc = Dataset(file, 'r')

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

lat2,lon2 = np.meshgrid(lat,lon)
time = nc.variables['time'][:]
z = nc.variables['bref'][:]



plt.figure(figsize=(14,14))

min_value = -5
max_value = 75
value_interval = 5
title_font_size = 22

cmin = min_value; cmax = max_value; cint = value_interval; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-103,-94]); ylim = np.array([32,37.5])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(); m.drawcountries()

cs = m.contourf(lon2,lat2,z[0,:,:].T, clevs, cmap = 'pyart_NWSRef', extend = 'both')

m.drawcounties()

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[dBZ]',size=title_font_size)
plt.title(r'KTLX 7/12 2351 UTC $Z_{H}$', size = title_font_size)
plt.show()

