#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 13:52:30 2022

@author: noahbrauer
"""



import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader



file = 'RAXPOL-20220928-234845-E6.0-Z.nc'

#Read in file

nc = Dataset(file, 'r')

#Extract attributes

azimuth = nc.variables['Azimuth'][:]
z = nc.variables['Intensity'][:]
gate = nc.variables['GateWidth'][:]
elevation = nc.variables['Elevation'][:]

z[z<0.1] = np.nan


max_range = 49800.003 - 30
gate_spacing = 30


radar_range = np.arange(0,max_range, step = gate_spacing)


x = radar_range*np.sin(np.deg2rad(azimuth))[:,None]
y = radar_range*np.cos(np.deg2rad(azimuth))[:,None]

RadarLongitude = -81.94765
RadarLatitude = 28.75541


#%%


proj = cartopy.crs.LambertConformal(central_longitude = RadarLongitude, central_latitude = RadarLatitude)
    
state_borders = cartopy.feature.NaturalEarthFeature(category='cultural', name = 'admin_1_states_provinces_lakes', scale = '50m', facecolor = 'none')

reader = shpreader.Reader('countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())    



fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = 0.; cmax = 60.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)



mesh = ax.pcolormesh(x,y,z, cmap = 'turbo', vmin = 0, vmax = 60)

ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='black')

distance_in_degrees = 0.3

ax.set_extent([RadarLongitude - distance_in_degrees, RadarLongitude + distance_in_degrees, RadarLatitude - distance_in_degrees, RadarLatitude + distance_in_degrees]) 



cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel('[dBZ]',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(18)
    
  
    
plt.title(r'RaXPol $6.0^{o}$ $Z_{H}$ 9/28 2348 UTC', size = 24) 





    


