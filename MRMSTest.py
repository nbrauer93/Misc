#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:58:17 2020

@author: noahbrauer
"""

import pygrib
import matplotlib.pyplot as plt
import numpy as np

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import pyart

grib = 'MRMS_SeamlessHSR_00.00_20200619-145600.grib2'

grbs = pygrib.open(grib)

#%%
grb = grbs.select()[0]
data = grb.values
#%%
# need to shift data grid longitudes from (0..360) to (-180..180)
lons = np.linspace(float(grb['longitudeOfFirstGridPointInDegrees']), float(grb['longitudeOfLastGridPointInDegrees']), int(grb['Ni']) )-360
lats = np.linspace(float(grb['latitudeOfFirstGridPointInDegrees']), float(grb['latitudeOfLastGridPointInDegrees']), int(grb['Nj']) )
grid_lon, grid_lat = np.meshgrid(lons, lats) #regularly spaced 2D grid



data[data<0]= np.nan

#print(np.unique(data[~np.isnan(data)]))

#%%

cmin = 0; cmax = 75; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(10,10))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-102.5,-94]); ylim = np.array([32,38.5])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(grid_lon,grid_lat,data,clevs,cmap='pyart_NWSRef',extend='both') 
m.drawcounties()
#m.pcolor(lon, lat, significant30, hatch='.', alpha = 0.)
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dBZ',size=28)

plt.title(r'MRMS $Z_{H}$ 6/19 1456 UTC ',name='Calibri',size=30)

x2star,y2star = m(-102.25,32.2)
plt.text(x2star,y2star, 'Plot by Noah Brauer', color = 'k')

plt.savefig('/Users/noahbrauer/Desktop/Radar_MRMS/Images/1456.png')
#plt.show(block=False) 