#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 17:11:08 2019

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob




import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

q_file = 'shum.mon.mean.nc'
u_file = 'uwnd.mon.mean.nc'
v_file = 'vwnd.mon.mean.nc'


ncq = Dataset(q_file, 'r')
ncu = Dataset(u_file, 'r')
ncv = Dataset(v_file, 'r')




narr = {}
lat = ncq.variables['lat'][:]
lon = ncq.variables['lon'][:] 

time = ncq.variables['time'][:]
timeUnits = ncq.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])

time_index_climo = np.where((narr['month']==9))[0]



narr2 = {}

time2 = ncu.variables['time'][:]
timeUnits2 = ncu.variables['time'].units
tmpDates2 = num2date(time2,timeUnits2,calendar='gregorian')
narr2['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates2])
narr2['day'] = np.asarray([d.day for d in narr2['date']])
narr2['month'] = np.asarray([d.month for d in narr2['date']])
narr2['year'] = np.asarray([d.year for d in narr2['date']])

time_index_climo_wind = np.where((narr2['month']==8))[0]





#Read in data

q = ncq.variables['shum'][time_index_climo, :,:,:]
uwnd = ncu.variables['uwnd'][time_index_climo_wind, :,:,:]
vwnd = ncv.variables['vwnd'][time_index_climo_wind, :,:,:]

pres = ncq.variables['level'][:]
g = 9.8

#Compute mean  with time (29,277,349) or height,lat,lon

q_mean = np.nanmean(q, axis = 0)

#Compute IVT

ivt_u = np.trapz(q_mean*uwnd,pres, axis = 1)*(-1/g)
ivt_v = np.trapz(q_mean*vwnd,pres, axis = 1)*(-1/g)

#Now compute mean IVT over time

ivt_u_mean = np.nanmean(ivt_u, axis = 0)
ivt_v_mean = np.nanmean(ivt_v, axis = 0)

def magnitude(u,v):
    mag = np.sqrt((u**2) +(v**2))
    return mag



ivt_mag_mean = magnitude(ivt_u_mean,ivt_v_mean)


#Compute standard deviations

ivt_u_std = np.nanstd(ivt_u, axis = 0)
ivt_v_std = np.nanstd(ivt_u, axis = 0)

ivt_mag_std = magnitude(ivt_u_std, ivt_v_std)


ivt_mag_nan = np.ones((277,349))*np.nan
ivt_u_nan = np.ones((277,349))*np.nan
ivt_v_nan = np.ones((277,349))*np.nan

for i in range(ivt_mag_nan.shape[0]):
    for j in range(ivt_mag_nan.shape[1]):
        
        if lon[i,j]>0:
            ivt_mag_nan[i,j] = np.nan
            ivt_u_nan[i,j] = np.nan
            ivt_v_nan[i,j] = np.nan
            
            
        else:
            ivt_mag_nan[i,j] = ivt_mag_mean[i,j]
            ivt_u_nan[i,j] = ivt_u_mean[i,j]
            ivt_v_nan[i,j] = ivt_v_mean[i,j]
            



#%%
#Now Plot


cmin = 0.; cmax = 4.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    
plt.figure(figsize=(10,6))
  

#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-130,-75]); ylim = np.array([15,50])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m.contourf(lon,lat,ivt_mag_nan,clevs,cmap='Greens',extend='both') 
widths = np.linspace(0, 2, lat.size)
#m.quiver(lon[::6,::6],lat[::6,::6],ivt_u[8,::6,::6],ivt_v[8,::6,::6], scale = 100)
m.quiver(lon[::4,::4],lat[::4,::4],ivt_u_nan[::4,::4],ivt_v_nan[::4,::4], scale = 150)


#m.drawcounties()

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel(r'kg $m^{-1}$ $s^{-1}$',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
   
plt.title('August Mean Integrated Vapor Transport (1979-2018 Climatology)', size = 20)


#x2star,y2star = m(-128,21)
#plt.text(x2star,y2star, 'Plot by Noah Brauer', weight = 'bold', color = 'b')


plt.show(block=False)











