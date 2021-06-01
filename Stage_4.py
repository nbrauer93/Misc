# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 12:14:38 2021

@author: noahb
"""

import numpy as np
import matplotlib.pyplot as plt


from netCDF4 import Dataset, num2date, MFDataset
import netCDF4 as nc


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap
import pyart


file = 'conus_20160417_24h.nc'


nc = Dataset(file, 'r')

#Import the radar attributes


precip  = nc.variables['A_PCP_GDS5_SFC_acc24h'][:]
lat = nc.variables['g5_lat_0'][:]
lon = nc.variables['g5_lon_1'][:]


precip[precip<=0] = np.nan

#%%

#Now let's plot


cmin = 0; cmax = 250; cint = 25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(10,10))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-102.5,-93]); ylim = np.array([28,32])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon,lat,precip,clevs,cmap='pyart_NWSRef',extend='both') 
m.drawcounties()


m.readshapefile("USCounties",'Watersheds2')

#m.pcolor(lon, lat, significant30, hatch='.', alpha = 0.)
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[mm]',size=28)

plt.title(r'Stage IV Rainfall 4/16 12Z - 4/17 12Z',name='Calibri',size=30)




