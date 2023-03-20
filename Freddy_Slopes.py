 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:42:26 2022

@author: noahbrauer
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj
from mpl_toolkits.basemap import Basemap
from shapely.geometry import LineString
from scipy import stats
from matplotlib.patches import Circle



#From Typhoon

#file = '2A.GPM.DPR.V8-20180723.20180911-S034010-E051243.025770.V06A.HDF5'


file = '2A.GPM.DPR.V9-20211125.20230208-S133017-E150250.050838.V07A.HDF5'


DPR = h5py.File(file, 'r')

lat = DPR['FS']['Latitude'][:,:]
lon = DPR['FS']['Longitude'][:,:]
z = DPR['FS']['SLV']['zFactorFinalESurface'][:] #nscan x nray (7934,49)
freezing = DPR['FS']['VER']['heightZeroDeg'][:]
dsd = DPR['FS']['SLV']['paramDSD'][:] #nscan x nray x nbin x DSDmoment (7934,49,176,2)

nw = dsd[:,:,:,0]/10
dm = dsd[:,:,:,1]

#heightStormTop

storm_top = DPR['FS']['PRE']['heightStormTop'][:]/1000

ku = DPR['FS']['SLV']['zFactorFinal'][:,:,:,0]
ku[ku<=12] = np.nan


storm_lat = -15.75
storm_lon =112.25

#Now add in shear vector

wind_dir = 31.7  #In degrees

#shear_dir = np.deg2rad(wind_dir) #Convert from degrees to radians


#Let's create a function to determine the end point of the arrow based on shear direction:

def arrow_endpoint(wind_dir):

    if 0<wind_dir<90: #NE shear

        angle_deg = 90 - wind_dir
        angle = np.deg2rad(angle_deg)
        dx = -1*np.cos(angle)
        dy = -1*np.sin(angle)

    if 90<wind_dir<180: #SE shear

        angle_deg = 180 - wind_dir
        angle = np.deg2rad(angle_deg)
        dx = -1*np.cos(angle)
        dy = np.sin(angle)

    if 180<wind_dir<270: #SW shear

        angle_deg = 270-wind_dir
        angle = np.deg2rad(angle_deg)
        dx = np.cos(angle)
        dy = np.sin(angle)

    if 270<wind_dir<360: #NW shear

        angle_deg = wind_dir - 270
        angle = np.deg2rad(angle_deg)
        dx = np.cos(angle)
        dy = -1*np.sin(angle)

    if wind_dir  == 0 or wind_dir == 360:

        dx = 0
        dy = -1

    if wind_dir == 90:

        dx = -1
        dy = 0

    if wind_dir == 180:

        dx = 0
        dy = 1

    if wind_dir == 270:

        dx = 1
        dy = 0


    return [dx,dy]


dx = arrow_endpoint(wind_dir)[0]
dy = arrow_endpoint(wind_dir)[1]

#print(dx)
#print(dy)


#Now compute and plot a perpendicular line


if dy != 0 :

    perp_slope = (-1*dx/dy)
    intercept = storm_lat+(dx/dy)*storm_lon

else:

    perp_slope = 0
    intercept = storm_lon

#print(intercept)
#print(perp_slope)


#Now set up two arbitrary arrays for lat and lon based off storm center:
#Line length refers to desired length of arrow (detertmine by radius of maximum wind?)

line_length = 0.75

lats = np.array([storm_lat-line_length, storm_lat + line_length])
lons = np.array([storm_lon-line_length, storm_lon + line_length])


#So now to plot line perpendicular to arrow:
# x = lons
# y = (perp_slope)*lons + intercept


#Now extend the arrow down past the point
#Arrow length is simply sqrt(dx**2 + dy**2)

arrow_length = np.sqrt((dx**2)+ (dy**2))

coords_lat = np.array([storm_lat,storm_lat - dy])
coords_lon = np.array([storm_lon,storm_lon- dx])


#%%

#Reshape lat and lon:

I,J = lat.shape

lat_reshape = lat.reshape(I*J, order = 'F')
lon_reshape = lon.reshape(I*J,order = 'F')


#Now import the function to partition by distance from storm center

def distance(x1,y1,x2,y2):

    dist = np.sqrt(((x2-x1)**2)+(y2-y1)**2) #Returns distance in degrees
    dist_km = dist*111

    return dist, dist_km

#Now compute distance from storm center:

distance_from_center = distance(storm_lon,storm_lat,lon_reshape,lat_reshape)[1]




#ur_quad = np.where(storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dy)&(storm_lon<=lon)&(lon<=storm_lon+np.abs(dx)))[0]

#%%
#Now partition KuPR (vertical profiles into shear-relative-quadrants)

#Ku has shape of latitude, angle (ray), height
#Let's contrain the latitude dimension to each shear-relative quadrant
#Index for each quadrant.
'''
from geopy.distance import geodesic
geodesic(kilometers=500).destination((start_lat,start_lon), 115)
'''


eyewall_dist = 100
inner_core_range = np.array([100,250])



def shear_quadrants(lon, lat, storm_lon, storm_lat, shear_dir, distance_from_center):

    if 0<shear_dir<90:

        dr_quad_eyewall = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon-np.abs(dx)-np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center<eyewall_dist))[0]
        dl_quad_eyewall = np.where((storm_lat-np.abs(dy) - np.abs(dx) <=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]
        ul_quad_eyewall = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon<=lon)&(lon<storm_lon+np.abs(dx) + np.abs (dy))&(distance_from_center<eyewall_dist))[0]
        ur_quad_eyewall = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+np.abs(dx))&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]

        dr_quad_inner = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon-np.abs(dx)-np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        dl_quad_inner = np.where((storm_lat-np.abs(dy) - np.abs(dx) <=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        ul_quad_inner = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon<=lon)&(lon<storm_lon+np.abs(dx) + np.abs (dy))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        ur_quad_inner = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+np.abs(dx))&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>45)&(distance_from_center<=250))[0]


    if 90<shear_dir<180:

        dr_quad_eyewall = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+np.abs(dx))&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]
        dl_quad_eyewall = np.where((storm_lat-np.abs(dx)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon-np.abs(dx)- np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center<eyewall_dist))[0]
        ul_quad_eyewall = np.where((storm_lat-np.abs(dy)- np.abs(dx)<=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dy)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]
        ur_quad_eyewall = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dx)+ np.abs(dx))&(storm_lon<=lon)&(lon<=storm_lon+np.abs(dx) + np.abs(dy))&(distance_from_center<eyewall_dist))[0]

        dr_quad_inner = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+np.abs(dx))&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        dl_quad_inner = np.where((storm_lat-np.abs(dx)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon-np.abs(dx)- np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        ul_quad_inner = np.where((storm_lat-np.abs(dy)- np.abs(dx)<=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dy)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        ur_quad_inner = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dx)+ np.abs(dx))&(storm_lon<=lon)&(lon<=storm_lon+np.abs(dx) + np.abs(dy))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]


    if 180<shear_dir<270:
        #Good
        dr_quad_eyewall = np.where((storm_lat-np.abs(dx)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon<=lon)&(lon<=storm_lon+np.abs(dx)+np.abs(dy))&(distance_from_center<eyewall_dist))[0]
        dl_quad_eyewall = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+np.abs(dx))&(storm_lon-np.abs(dx)-np.abs(dy)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]
        ul_quad_eyewall = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dy)+ np.abs(dx))&(storm_lon-np.abs(dx)- np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center<eyewall_dist))[0]
        ur_quad_eyewall = np.where((storm_lat-np.abs(dy) - np.abs(dx)<=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]

        dr_quad_inner = np.where((storm_lat-np.abs(dx)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon<=lon)&(lon<=storm_lon+np.abs(dx)+np.abs(dy))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        dl_quad_inner = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+np.abs(dx))&(storm_lon-np.abs(dx)-np.abs(dy)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        ul_quad_inner = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dy)+ np.abs(dx))&(storm_lon-np.abs(dx)- np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        ur_quad_inner = np.where((storm_lat-np.abs(dy) - np.abs(dx)<=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]


    if 270<shear_dir<360:
        #Good
        dr_quad_eyewall = np.where((storm_lat-np.abs(dy)- np.abs(dx)<=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]
        dl_quad_eyewall = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dx))&(storm_lon<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]
        ul_quad_eyewall = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+ np.abs(dx))&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center<eyewall_dist))[0]
        ur_quad_eyewall = np.where((storm_lat-np.abs(dx)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon-np.abs(dx)- np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center<eyewall_dist))[0]

        dr_quad_inner = np.where((storm_lat-np.abs(dy)- np.abs(dx)<=lat)&(lat<=storm_lat)&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        dl_quad_inner = np.where((storm_lat-np.abs(dy)<=lat)&(lat<=storm_lat+np.abs(dx))&(storm_lon<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]
        ul_quad_inner = np.where((storm_lat<=lat)&(lat<=storm_lat+np.abs(dy)+ np.abs(dx))&(storm_lon-np.abs(dx)<=lon)&(lon<=storm_lon+np.abs(dx))&(distance_from_center>inner_core_range[0])&(distance_from_center<=[1]))[0]
        ur_quad_inner = np.where((storm_lat-np.abs(dx)<=lat)&(lat<=storm_lat+np.abs(dy))&(storm_lon-np.abs(dx)- np.abs(dy)<=lon)&(lon<=storm_lon)&(distance_from_center>inner_core_range[0])&(distance_from_center<=inner_core_range[1]))[0]


    return dr_quad_eyewall,dl_quad_eyewall,ul_quad_eyewall,ur_quad_eyewall,dr_quad_inner,dl_quad_inner,ul_quad_inner,ur_quad_inner

#%%

#Now index ku into each quadrant:


dr_index_eyewall = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[0]
dl_index_eyewall = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[1]
ul_index_eyewall = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[2]
ur_index_eyewall = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[3]

'''
dr_index_inner = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[4]
dl_index_inner = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[5]
ul_index_inner = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[6]
ur_index_inner = shear_quadrants(lon_reshape,lat_reshape,storm_lon,storm_lat,wind_dir,distance_from_center)[7]
'''





#Designate a Z-direction
x_bins = np.arange(10,60,1)
y_bins = np.arange(0,21.988032,0.124932)






#plt.plot(ku_1.T,y_bins[::-1])
#plt.show()

#ku_dr_quad has shape (461,49,176) or footprint, angle, height








#Reshape KuPR to be NINDe x height


I,J,K = ku.shape

ku_reshape = ku.reshape(I*J,K, order = 'F')
freezing_reshape = freezing.reshape(I*J,order = 'F')/1000 # Convert to meters
dm_reshape  = dm.reshape(I*J,K ,order = 'F')
nw_reshape = nw.reshape(I*J,K, order = 'F')
top_reshape = storm_top.reshape(I*J, order = 'F')

top_reshape[top_reshape<0] = np.nan


#Normalize w.r.t to the zero degree isotherm
#Altitude - 0deg isotherm
#Loop through each freezing level height for all footprints and all angles

freezing_norm = []


for i in range(len(freezing_reshape)):

    freezing_normalized = y_bins - freezing_reshape[i]
    freezing_norm.append(freezing_normalized)


freezing_norm = np.asarray(freezing_norm)



#Now index for each shear relative quadrant and annulus:


ku_dl_eyewall = ku_reshape[dl_index_eyewall,::-1]
ku_dr_eyewall = ku_reshape[dr_index_eyewall,::-1]
ku_ul_eyewall = ku_reshape[ul_index_eyewall,::-1]
ku_ur_eyewall = ku_reshape[ur_index_eyewall,::-1]
'''
ku_dl_inner = ku_reshape[dl_index_inner,::-1]
ku_dr_inner = ku_reshape[dr_index_inner,::-1]
ku_ul_inner = ku_reshape[ul_index_inner,::-1]
ku_ur_inner = ku_reshape[ur_index_inner,::-1]


nw_dl_eyewall = nw_reshape[dl_index_eyewall,::-1]
nw_dr_eyewall = nw_reshape[dr_index_eyewall,::-1]
nw_ul_eyewall = nw_reshape[ul_index_eyewall,::-1]
nw_ur_eyewall = nw_reshape[ur_index_eyewall,::-1]

nw_dl_inner = nw_reshape[dl_index_inner,::-1]
nw_dr_inner = nw_reshape[dr_index_inner,::-1]
nw_ul_inner = nw_reshape[ul_index_inner,::-1]
nw_ur_inner = nw_reshape[ur_index_inner,::-1]


dm_dl_eyewall = dm_reshape[dl_index_eyewall,::-1]
dm_dr_eyewall = dm_reshape[dr_index_eyewall,::-1]
dm_ul_eyewall = dm_reshape[ul_index_eyewall,::-1]
dm_ur_eyewall = dm_reshape[ur_index_eyewall,::-1]

dm_dl_inner = dm_reshape[dl_index_inner,::-1]
dm_dr_inner = dm_reshape[dr_index_inner,::-1]
dm_ul_inner = dm_reshape[ul_index_inner,::-1]
dm_ur_inner = dm_reshape[ur_index_inner,::-1]
'''

top_dl_eyewall = top_reshape[dl_index_eyewall]
top_dr_eyewall = top_reshape[dr_index_eyewall]
top_ul_eyewall = top_reshape[ul_index_eyewall]
top_ur_eyewall = top_reshape[ur_index_eyewall]
'''
top_dl_inner = top_reshape[dl_index_inner]
top_dr_inner = top_reshape[dr_index_inner]
top_ul_inner = top_reshape[ul_index_inner]
top_ur_inner = top_reshape[ur_index_inner]
'''




freezing_dl_eyewall = freezing_norm[dl_index_eyewall]
freezing_dr_eyewall = freezing_norm[dr_index_eyewall]
freezing_ul_eyewall = freezing_norm[ul_index_eyewall]
freezing_ur_eyewall = freezing_norm[ur_index_eyewall]
'''
freezing_dl_inner = freezing_norm[dl_index_inner]
freezing_dr_inner = freezing_norm[dr_index_inner]
freezing_ul_inner = freezing_norm[ul_index_inner]
freezing_ur_inner = freezing_norm[ur_index_inner]

'''


#Now we need to compute slope both within the warm cloud layer and the ice-phase layer
#Mask everything greater than normalized altitude of -1 km abd everything less than -3 km

ku_dl_eyewall_warm = ku_dl_eyewall.copy()
ku_dl_eyewall_warm[(freezing_dl_eyewall>-1) | (freezing_dl_eyewall<-3)] = -9999
ku_dl_eyewall_warm = np.ma.masked_where(ku_dl_eyewall_warm<0,ku_dl_eyewall_warm)

ku_dr_eyewall_warm = ku_dr_eyewall.copy()
ku_dr_eyewall_warm[(freezing_dr_eyewall>-1) | (freezing_dr_eyewall<-3)] = -9999
ku_dr_eyewall_warm = np.ma.masked_where(ku_dr_eyewall_warm<0,ku_dr_eyewall_warm)

ku_ul_eyewall_warm = ku_ul_eyewall.copy()
ku_ul_eyewall_warm[(freezing_ul_eyewall>-1) | (freezing_ul_eyewall<-3)] = -9999
ku_ul_eyewall_warm = np.ma.masked_where(ku_ul_eyewall_warm<0,ku_ul_eyewall_warm)

ku_ur_eyewall_warm = ku_ur_eyewall.copy()
ku_ur_eyewall_warm[(freezing_ur_eyewall>-1) | (freezing_ur_eyewall<-3)] = -9999
ku_ur_eyewall_warm = np.ma.masked_where(ku_ur_eyewall_warm<0,ku_ur_eyewall_warm)

'''
ku_dl_inner_warm = ku_dl_inner.copy()
ku_dl_inner_warm[(freezing_dl_inner>-1) | (freezing_dl_inner<-3)] = -9999
ku_dl_inner_warm = np.ma.masked_where(ku_dl_inner_warm<0,ku_dl_inner_warm)

ku_dr_inner_warm = ku_dr_inner.copy()
ku_dr_inner_warm[(freezing_dr_inner>-1) | (freezing_dr_inner<-3)] = -9999
ku_dr_inner_warm = np.ma.masked_where(ku_dr_inner_warm<0,ku_dr_inner_warm)

ku_ul_inner_warm = ku_ul_inner.copy()
ku_ul_inner_warm[(freezing_ul_inner>-1) | (freezing_ul_inner<-3)] = -9999
ku_ul_inner_warm = np.ma.masked_where(ku_ul_inner_warm<0,ku_ul_inner_warm)

ku_ur_inner_warm = ku_ur_inner.copy()
ku_ur_inner_warm[(freezing_ur_inner>-1) | (freezing_ur_inner<-3)] = -9999
ku_ur_inner_warm = np.ma.masked_where(ku_ur_inner_warm<0,ku_ur_inner_warm)



#Now for Nw and dm


dm_dl_eyewall_warm = dm_dl_eyewall.copy()
dm_dl_eyewall_warm[(freezing_dl_eyewall>-1) | (freezing_dl_eyewall<-3)] = -9999
dm_dl_eyewall_warm = np.ma.masked_where(dm_dl_eyewall_warm<0,dm_dl_eyewall_warm)

dm_dr_eyewall_warm = dm_dr_eyewall.copy()
dm_dr_eyewall_warm[(freezing_dr_eyewall>-1) | (freezing_dr_eyewall<-3)] = -9999
dm_dr_eyewall_warm = np.ma.masked_where(dm_dr_eyewall_warm<0,dm_dr_eyewall_warm)

dm_ul_eyewall_warm = dm_ul_eyewall.copy()
dm_ul_eyewall_warm[(freezing_ul_eyewall>-1) | (freezing_ul_eyewall<-3)] = -9999
dm_ul_eyewall_warm = np.ma.masked_where(dm_ul_eyewall_warm<0,dm_ul_eyewall_warm)

dm_ur_eyewall_warm = dm_ur_eyewall.copy()
dm_ur_eyewall_warm[(freezing_ur_eyewall>-1) | (freezing_ur_eyewall<-3)] = -9999
dm_ur_eyewall_warm = np.ma.masked_where(dm_ur_eyewall_warm<0,dm_ur_eyewall_warm)


dm_dl_inner_warm = dm_dl_inner.copy()
dm_dl_inner_warm[(freezing_dl_inner>-1) | (freezing_dl_inner<-3)] = -9999
dm_dl_inner_warm = np.ma.masked_where(dm_dl_inner_warm<0,dm_dl_inner_warm)

dm_dr_inner_warm = dm_dr_inner.copy()
dm_dr_inner_warm[(freezing_dr_inner>-1) | (freezing_dr_inner<-3)] = -9999
dm_dr_inner_warm = np.ma.masked_where(dm_dr_inner_warm<0,dm_dr_inner_warm)

dm_ul_inner_warm = dm_ul_inner.copy()
dm_ul_inner_warm[(freezing_ul_inner>-1) | (freezing_ul_inner<-3)] = -9999
dm_ul_inner_warm = np.ma.masked_where(dm_ul_inner_warm<0,dm_ul_inner_warm)

dm_ur_inner_warm = dm_ur_inner.copy()
dm_ur_inner_warm[(freezing_ur_inner>-1) | (freezing_ur_inner<-3)] = -9999
dm_ur_inner_warm = np.ma.masked_where(dm_ur_inner_warm<0,dm_ur_inner_warm)





nw_dl_eyewall_warm = nw_dl_eyewall.copy()
nw_dl_eyewall_warm[(freezing_dl_eyewall>-1) | (freezing_dl_eyewall<-3)] = -9999
nw_dl_eyewall_warm = np.ma.masked_where(nw_dl_eyewall_warm<0,nw_dl_eyewall_warm)

nw_dr_eyewall_warm = nw_dr_eyewall.copy()
nw_dr_eyewall_warm[(freezing_dr_eyewall>-1) | (freezing_dr_eyewall<-3)] = -9999
nw_dr_eyewall_warm = np.ma.masked_where(nw_dr_eyewall_warm<0,nw_dr_eyewall_warm)

nw_ul_eyewall_warm = nw_ul_eyewall.copy()
nw_ul_eyewall_warm[(freezing_ul_eyewall>-1) | (freezing_ul_eyewall<-3)] = -9999
nw_ul_eyewall_warm = np.ma.masked_where(nw_ul_eyewall_warm<0,nw_ul_eyewall_warm)

nw_ur_eyewall_warm = nw_ur_eyewall.copy()
nw_ur_eyewall_warm[(freezing_ur_eyewall>-1) | (freezing_ur_eyewall<-3)] = -9999
nw_ur_eyewall_warm = np.ma.masked_where(nw_ur_eyewall_warm<0,nw_ur_eyewall_warm)


nw_dl_inner_warm = nw_dl_inner.copy()
nw_dl_inner_warm[(freezing_dl_inner>-1) | (freezing_dl_inner<-3)] = -9999
nw_dl_inner_warm = np.ma.masked_where(nw_dl_inner_warm<0,nw_dl_inner_warm)

nw_dr_inner_warm = nw_dr_inner.copy()
nw_dr_inner_warm[(freezing_dr_inner>-1) | (freezing_dr_inner<-3)] = -9999
nw_dr_inner_warm = np.ma.masked_where(nw_dr_inner_warm<0,nw_dr_inner_warm)

nw_ul_inner_warm = nw_ul_inner.copy()
nw_ul_inner_warm[(freezing_ul_inner>-1) | (freezing_ul_inner<-3)] = -9999
nw_ul_inner_warm = np.ma.masked_where(nw_ul_inner_warm<0,nw_ul_inner_warm)

nw_ur_inner_warm = nw_ur_inner.copy()
nw_ur_inner_warm[(freezing_ur_inner>-1) | (freezing_ur_inner<-3)] = -9999
nw_ur_inner_warm = np.ma.masked_where(nw_ur_inner_warm<0,nw_ur_inner_warm)

'''


#Now do the linear regression within the warm cloud layer to compute the slope




slope_warm_dl_eyewall = []


for i in range(ku_dl_eyewall.shape[0]):

    slope_dl_eyewall = stats.mstats.linregress(ku_dl_eyewall_warm[i,:],freezing_dl_eyewall[i,:])[0]
    slope_warm_dl_eyewall.append(slope_dl_eyewall)




slope_warm_dl_eyewall = np.asarray(slope_warm_dl_eyewall)




slope_warm_dr_eyewall = []


for i in range(ku_dr_eyewall.shape[0]):

    slope_dr_eyewall = stats.mstats.linregress(ku_dr_eyewall_warm[i,:],freezing_dr_eyewall[i,:])[0]
    slope_warm_dr_eyewall.append(slope_dr_eyewall)




slope_warm_dr_eyewall = np.asarray(slope_warm_dr_eyewall)




slope_warm_ul_eyewall = []


for i in range(ku_ul_eyewall.shape[0]):

    slope_ul_eyewall = stats.mstats.linregress(ku_ul_eyewall_warm[i,:],freezing_ul_eyewall[i,:])[0]
    slope_warm_ul_eyewall.append(slope_ul_eyewall)




slope_warm_ul_eyewall = np.asarray(slope_warm_ul_eyewall)




slope_warm_ur_eyewall = []


for i in range(ku_ur_eyewall.shape[0]):

    slope_ur_eyewall = stats.mstats.linregress(ku_ur_eyewall_warm[i,:],freezing_ur_eyewall[i,:])[0]
    slope_warm_ur_eyewall.append(slope_ur_eyewall)




slope_warm_ur_eyewall = np.asarray(slope_warm_ur_eyewall)




#Now for the inner core


'''
slope_warm_dl_inner = []
slope_nw_dl_inner = []
slope_dm_dl_inner = []

for i in range(ku_dl_inner.shape[0]):

    slope_dl_inner = stats.mstats.linregress(ku_dl_inner_warm[i,:],freezing_dl_inner[i,:])[0]
    slope_warm_dl_inner.append(slope_dl_inner)

    slope_dl_inner_nw = stats.mstats.linregress(nw_dl_inner_warm[i,:],freezing_dl_inner[i,:])[0]
    slope_nw_dl_inner.append(slope_dl_inner_nw)

    slope_dl_inner_dm = stats.mstats.linregress(dm_dl_inner_warm[i,:],freezing_dl_inner[i,:])[0]
    slope_dm_dl_inner.append(slope_dl_inner_dm)


slope_warm_dl_inner = np.asarray(slope_warm_dl_inner)
slope_nw_dl_inner = np.asarray(slope_nw_dl_inner)
slope_dm_dl_inner = np.asarray(slope_dm_dl_inner)





slope_warm_dr_inner = []
slope_nw_dr_inner = []
slope_dm_dr_inner = []

for i in range(ku_dr_inner.shape[0]):

    slope_dr_inner = stats.mstats.linregress(ku_dr_inner_warm[i,:],freezing_dr_inner[i,:])[0]
    slope_warm_dr_inner.append(slope_dr_inner)

    slope_dr_inner_nw = stats.mstats.linregress(nw_dr_inner_warm[i,:],freezing_dr_inner[i,:])[0]
    slope_nw_dr_inner.append(slope_dr_inner_nw)


    slope_dr_inner_dm = stats.mstats.linregress(dm_dr_inner_warm[i,:],freezing_dr_inner[i,:])[0]
    slope_dm_dr_inner.append(slope_dr_inner_dm)


slope_warm_dl_inner = np.asarray(slope_warm_dl_inner)
slope_nw_dl_inner = np.asarray(slope_nw_dl_inner)
slope_dm_dl_inner = np.asarray(slope_dm_dl_inner)



slope_warm_ul_inner = []
slope_nw_ul_inner = []
slope_dm_ul_inner = []

for i in range(ku_ul_inner.shape[0]):

    slope_ul_inner = stats.mstats.linregress(ku_ul_inner_warm[i,:],freezing_ul_inner[i,:])[0]
    slope_warm_ul_inner.append(slope_ul_inner)

    slope_ul_inner_nw = stats.mstats.linregress(nw_ul_inner_warm[i,:],freezing_ul_inner[i,:])[0]
    slope_nw_ul_inner.append(slope_ul_inner_nw)


    slope_ul_inner_dm = stats.mstats.linregress(dm_ul_inner_warm[i,:],freezing_ul_inner[i,:])[0]
    slope_dm_ul_inner.append(slope_ul_inner_dm)


slope_warm_ul_inner = np.asarray(slope_warm_ul_inner)
slope_nw_ul_inner = np.asarray(slope_nw_ul_inner)
slope_dm_ul_inner = np.asarray(slope_dm_ul_inner)



slope_warm_ur_inner = []
slope_nw_ur_inner = []
slope_dm_ur_inner = []

for i in range(ku_ur_inner.shape[0]):

    slope_ur_inner = stats.mstats.linregress(ku_ur_inner_warm[i,:],freezing_ur_inner[i,:])[0]
    slope_warm_ur_inner.append(slope_ur_inner)

    slope_ur_inner_nw = stats.mstats.linregress(nw_ur_inner_warm[i,:],freezing_ur_inner[i,:])[0]
    slope_nw_ur_inner.append(slope_ur_inner_nw)

    slope_ur_inner_dm = stats.mstats.linregress(dm_ur_inner_warm[i,:],freezing_ur_inner[i,:])[0]
    slope_dm_ur_inner.append(slope_ur_inner_dm)


slope_warm_ur_inner = np.asarray(slope_warm_ur_inner)
slope_nw_ur_inner = np.asarray(slope_nw_ur_inner)
slope_dm_ur_inner = np.asarray(slope_dm_ur_inner)

'''

#Let's plot PDFs for all



'''
#Now for Nw and Dm


sns.distplot(slope_nw_dl_eyewall, hist = False, kde = True, bins = 50, color = 'darkblue', label = 'Downshear Left')
sns.distplot(slope_nw_dr_eyewall, hist = False, kde = True, bins = 50, color = 'red', label = 'Downshear Right')
sns.distplot(slope_nw_ul_eyewall, hist = False, kde = True, bins = 50, color = 'k', label = 'Upshear Left')
sns.distplot(slope_nw_ur_eyewall, hist = False, kde = True, bins = 50, color = 'green', label = 'Upshear Right')
plt.legend()
plt.title(r' $log_{10}(N_{W})$ Slope in Liquid Phase in Eyewall', size = 16)
plt.xlabel(r'[$m^{-3}$ $mm^{-1}$ $km^{-1}$]', size = 16)
plt.ylabel('Density',size = 16)

plt.show()





sns.distplot(slope_nw_dl_inner, hist = False, kde = True, bins = 50, color = 'darkblue', label = 'Downshear Left')
sns.distplot(slope_nw_dr_inner, hist = False, kde = True, bins = 50, color = 'red', label = 'Downshear Right')
sns.distplot(slope_nw_ul_inner, hist = False, kde = True, bins = 50, color = 'k', label = 'Upshear Left')
sns.distplot(slope_nw_ur_inner, hist = False, kde = True, bins = 50, color = 'green', label = 'Upshear Right')
plt.legend()
plt.title('$log_{10}(N_{W})$ Slope in Liquid Phase in Inner Core', size = 16)
plt.xlabel(r'[$m^{-3}$ $mm^{-1}$ $km^{-1}$]', size = 16)
plt.ylabel('Density',size = 16)
plt.show()



sns.distplot(slope_dm_dl_eyewall, hist = False, kde = True, bins = 50, color = 'darkblue', label = 'Downshear Left')
sns.distplot(slope_dm_dr_eyewall, hist = False, kde = True, bins = 50, color = 'red', label = 'Downshear Right')
sns.distplot(slope_dm_ul_eyewall, hist = False, kde = True, bins = 50, color = 'k', label = 'Upshear Left')
sns.distplot(slope_dm_ur_eyewall, hist = False, kde = True, bins = 50, color = 'green', label = 'Upshear Right')
plt.legend()
plt.title(r' $D_{M}$ Slope in Liquid Phase in Eyewall', size = 16)
plt.xlabel(r'[mm/km]', size = 16)
plt.ylabel('Density',size = 16)

plt.show()





sns.distplot(slope_dm_dl_inner, hist = False, kde = True, bins = 50, color = 'darkblue', label = 'Downshear Left')
sns.distplot(slope_dm_dr_inner, hist = False, kde = True, bins = 50, color = 'red', label = 'Downshear Right')
sns.distplot(slope_dm_ul_inner, hist = False, kde = True, bins = 50, color = 'k', label = 'Upshear Left')
sns.distplot(slope_dm_ur_inner, hist = False, kde = True, bins = 50, color = 'green', label = 'Upshear Right')
plt.legend()
plt.title('$D_{M}$ Slope in Liquid Phase in Inner Core', size = 16)
plt.xlabel(r'[mm/km]', size = 16)
plt.ylabel('Density',size = 16)
plt.show()
'''
#%%
#Lastly for echo top height









#%%

#Now let's do this for the ice-phase


#Normalized altitude of 1 to 4 km W.R.T freezing level

ku_dl_eyewall_ice = ku_dl_eyewall.copy()
ku_dl_eyewall_ice[(freezing_dl_eyewall>4) | (freezing_dl_eyewall<1)] = -9999
ku_dl_eyewall_ice = np.ma.masked_where(ku_dl_eyewall_ice<0,ku_dl_eyewall_ice)

ku_dr_eyewall_ice = ku_dr_eyewall.copy()
ku_dr_eyewall_ice[(freezing_dr_eyewall>4) | (freezing_dr_eyewall<1)] = -9999
ku_dr_eyewall_ice = np.ma.masked_where(ku_dr_eyewall_ice<0,ku_dr_eyewall_ice)

ku_ul_eyewall_ice = ku_ul_eyewall.copy()
ku_ul_eyewall_ice[(freezing_ul_eyewall>4) | (freezing_ul_eyewall<1)] = -9999
ku_ul_eyewall_ice = np.ma.masked_where(ku_ul_eyewall_ice<0,ku_ul_eyewall_ice)

ku_ur_eyewall_ice = ku_ur_eyewall.copy()
ku_ur_eyewall_ice[(freezing_ur_eyewall>4) | (freezing_ur_eyewall<1)] = -9999
ku_ur_eyewall_ice = np.ma.masked_where(ku_ur_eyewall_ice<0,ku_ur_eyewall_ice)
'''

ku_dl_inner_ice = ku_dl_inner.copy()
ku_dl_inner_ice[(freezing_dl_inner>4) | (freezing_dl_inner<1)] = -9999
ku_dl_inner_ice = np.ma.masked_where(ku_dl_inner_ice<0,ku_dl_inner_ice)

ku_dr_inner_ice = ku_dr_inner.copy()
ku_dr_inner_ice[(freezing_dr_inner>4) | (freezing_dr_inner<1)] = -9999
ku_dr_inner_ice = np.ma.masked_where(ku_dr_inner_ice<0,ku_dr_inner_ice)

ku_ul_inner_ice = ku_ul_inner.copy()
ku_ul_inner_ice[(freezing_ul_inner>4) | (freezing_ul_inner<1)] = -9999
ku_ul_inner_ice = np.ma.masked_where(ku_ul_inner_ice<0,ku_ul_inner_ice)

ku_ur_inner_ice = ku_ur_inner.copy()
ku_ur_inner_ice[(freezing_ur_inner>4) | (freezing_ur_inner<1)] = -9999
ku_ur_inner_ice = np.ma.masked_where(ku_ur_inner_ice<0,ku_ur_inner_ice)
'''

#%%

#Now perform the linear regression


slope_ice_dl_eyewall = []

for i in range(ku_dl_eyewall.shape[0]):

    slope_dl_eyewall = stats.mstats.linregress(ku_dl_eyewall_ice[i,:],freezing_dl_eyewall[i,:])[0]
    slope_ice_dl_eyewall.append(slope_dl_eyewall)



slope_ice_dr_eyewall = []

for i in range(ku_dr_eyewall.shape[0]):

    slope_dr_eyewall = stats.mstats.linregress(ku_dr_eyewall_ice[i,:],freezing_dr_eyewall[i,:])[0]
    slope_ice_dr_eyewall.append(slope_dr_eyewall)




slope_ice_ul_eyewall = []

for i in range(ku_ul_eyewall.shape[0]):

    slope_ul_eyewall = stats.mstats.linregress(ku_ul_eyewall_ice[i,:],freezing_ul_eyewall[i,:])[0]
    slope_ice_ul_eyewall.append(slope_ul_eyewall)





slope_ice_ur_eyewall = []

for i in range(ku_ur_eyewall.shape[0]):

    slope_ur_eyewall = stats.mstats.linregress(ku_ur_eyewall_ice[i,:],freezing_ur_eyewall[i,:])[0]
    slope_ice_ur_eyewall.append(slope_ur_eyewall)



#Now for the inner core

'''

slope_ice_dl_inner = []

for i in range(ku_dl_inner.shape[0]):

    slope_dl_inner = stats.mstats.linregress(ku_dl_inner_ice[i,:],freezing_dl_inner[i,:])[0]
    slope_ice_dl_inner.append(slope_dl_inner)



slope_ice_dr_inner = []

for i in range(ku_dr_inner.shape[0]):

    slope_dr_inner = stats.mstats.linregress(ku_dr_inner_ice[i,:],freezing_dr_inner[i,:])[0]
    slope_ice_dr_inner.append(slope_dr_inner)




slope_ice_ul_inner = []

for i in range(ku_ul_inner.shape[0]):

    slope_ul_inner = stats.mstats.linregress(ku_ul_inner_ice[i,:],freezing_ul_inner[i,:])[0]
    slope_ice_ul_inner.append(slope_ul_inner)





slope_ice_ur_inner = []

for i in range(ku_ur_inner.shape[0]):

    slope_ur_inner = stats.mstats.linregress(ku_ur_inner_ice[i,:],freezing_ur_inner[i,:])[0]
    slope_ice_ur_inner.append(slope_ur_inner)
#%%
'''



#%%


#Output all files here


slope_warm_dl_eyewall_file = np.save('0208_slope_warm_dl_eyewall.npy', slope_warm_dl_eyewall, allow_pickle = True)
slope_warm_dr_eyewall_file = np.save('0208_slope_warm_dr_eyewall.npy', slope_warm_dr_eyewall, allow_pickle = True)
slope_warm_ur_eyewall_file = np.save('0208_slope_warm_ur_eyewall.npy', slope_warm_ur_eyewall, allow_pickle = True)
slope_warm_ul_eyewall_file = np.save('0208_slope_warm_ul_eyewall.npy', slope_warm_ul_eyewall, allow_pickle = True)

'''
slope_nw_dl_eyewall_file = np.save('1_slope_nw_dl_eyewall_NIND.npy', slope_nw_dl_eyewall, allow_pickle = True)
slope_nw_dr_eyewall_file = np.save('1_slope_nw_dr_eyewall_NIND.npy', slope_nw_dr_eyewall, allow_pickle = True)
slope_nw_ur_eyewall_file = np.save('1_slope_nw_ur_eyewall_NIND.npy', slope_nw_ur_eyewall, allow_pickle = True)
slope_nw_ul_eyewall_file = np.save('1_slope_nw_ul_eyewall_NIND.npy', slope_nw_ul_eyewall, allow_pickle = True)

slope_dm_dl_eyewall_file = np.save('1_slope_dm_dl_eyewall_NIND.npy', slope_dm_dl_eyewall, allow_pickle = True)
slope_dm_dr_eyewall_file = np.save('1_slope_dm_dr_eyewall_NIND.npy', slope_dm_dr_eyewall, allow_pickle = True)
slope_dm_ur_eyewall_file = np.save('1_slope_dm_ur_eyewall_NIND.npy', slope_dm_ur_eyewall, allow_pickle = True)
slope_dm_ul_eyewall_file = np.save('1_slope_dm_ul_eyewall_NIND.npy', slope_dm_ul_eyewall, allow_pickle = True)
'''

top_dl_eyewall_file =  np.save('0208_top_dl_eyewall.npy', top_dl_eyewall, allow_pickle = True)
top_dr_eyewall_file =  np.save('0208_top_dr_eyewall.npy', top_dr_eyewall, allow_pickle = True)
top_ul_eyewall_file =  np.save('0208_top_ul_eyewall.npy', top_ul_eyewall, allow_pickle = True)
top_ur_eyewall_file =  np.save('0208_top_ur_eyewall.npy', top_ur_eyewall, allow_pickle = True)


'''
slope_warm_dl_inner_file = np.save('1_slope_warm_dl_inner_NIND.npy', slope_warm_dl_inner, allow_pickle = True)
slope_warm_dr_inner_file = np.save('1_slope_warm_dr_inner_NIND.npy', slope_warm_dr_inner, allow_pickle = True)
slope_warm_ur_inner_file = np.save('1_slope_warm_ur_inner_NIND.npy', slope_warm_ur_inner, allow_pickle = True)
slope_warm_ul_inner_file = np.save('1_slope_warm_ul_inner_NIND.npy', slope_warm_ul_inner, allow_pickle = True)
'''
'''
slope_nw_dl_inner_file = np.save('1_slope_nw_dl_inner_NIND.npy', slope_nw_dl_inner, allow_pickle = True)
slope_nw_dr_inner_file = np.save('1_slope_nw_dr_inner_NIND.npy', slope_nw_dr_inner, allow_pickle = True)
slope_nw_ur_inner_file = np.save('1_slope_nw_ur_inner_NIND.npy', slope_nw_ur_inner, allow_pickle = True)
slope_nw_ul_inner_file = np.save('1_slope_nw_ul_inner_NIND.npy', slope_nw_ul_inner, allow_pickle = True)

slope_dm_dl_inner_file = np.save('1_slope_dm_dl_inner_NIND.npy', slope_dm_dl_inner, allow_pickle = True)
slope_dm_dr_inner_file = np.save('1_slope_dm_dr_inner_NIND.npy', slope_dm_dr_inner, allow_pickle = True)
slope_dm_ur_inner_file = np.save('1_slope_dm_ur_inner_NIND.npy', slope_dm_ur_inner, allow_pickle = True)
slope_dm_ul_inner_file = np.save('1_slope_dm_ul_inner_NIND.npy', slope_dm_ul_inner, allow_pickle = True)
'''
'''
top_dl_inner_file =  np.save('1_top_dl_inner_NIND.npy', top_dl_inner, allow_pickle = True)
top_dr_inner_file =  np.save('1_top_dr_inner_NIND.npy', top_dr_inner, allow_pickle = True)
top_ul_inner_file =  np.save('1_top_ul_inner_NIND.npy', top_ul_inner, allow_pickle = True)
top_ur_inner_file =  np.save('1_top_ur_inner_NIND.npy', top_ur_inner, allow_pickle = True)

'''

slope_ice_dl_eyewall_file = np.save('0208_slope_ice_dl_eyewall.npy', slope_ice_dl_eyewall, allow_pickle = True)
slope_ice_dr_eyewall_file = np.save('0208_slope_ice_dr_eyewall.npy', slope_ice_dr_eyewall, allow_pickle = True)
slope_ice_ur_eyewall_file = np.save('0208_slope_ice_ur_eyewall.npy', slope_ice_ur_eyewall, allow_pickle = True)
slope_ice_ul_eyewall_file = np.save('0208_slope_ice_ul_eyewall.npy', slope_ice_ul_eyewall, allow_pickle = True)

'''
slope_ice_dl_inner_file = np.save('1_slope_ice_dl_inner_NIND.npy', slope_ice_dl_inner, allow_pickle = True)
slope_ice_dr_inner_file = np.save('1_slope_ice_dr_inner_NIND.npy', slope_ice_dr_inner, allow_pickle = True)
slope_ice_ur_inner_file = np.save('1_slope_ice_ur_inner_NIND.npy', slope_ice_ur_inner, allow_pickle = True)
slope_ice_ul_inner_file = np.save('1_slope_ice_ul_inner_NIND.npy', slope_ice_ul_inner, allow_pickle = True)
'''

'''
#%%


#Loop through and plot all vertical profiles in a designated quadrant


#Now plot on a lat-lon grid

z[z<12] = np.nan


cmin = 12.; cmax = 70.; cint = 2; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)

plt.figure(figsize=(14,14))


xlim = np.array([storm_lon-4,storm_lon+4]); ylim = np.array([storm_lat-4,storm_lat+4])


label_size = 22


m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,z[:,:,0],clevs,cmap='turbo',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
parallels = np.arange(-90,90,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(0,360, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])

cbar = plt.colorbar(fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.clim(12,70)

x2star,y2star = m(storm_lon,storm_lat)
#m.plot(x2star,y2star,markersize=26, color = 'k')
m.scatter(x2star,y2star,s=200, color = 'k')


arrow_length_x = dx
arrow_length_y = dy


x,y = m(storm_lon,storm_lat)
x2,y2 = m(dx,dy)
plt.arrow(x,y,x2,y2, width = 0.008, head_width = 0.2, head_length = 0.2, color = 'k')

x3,y3 = lons, (perp_slope)*lons + intercept
m.plot(x3,y3,linewidth = 2, color = 'k')

#Now draw another line extending down from the storm_lon, storm_lat point

x4,y4 = coords_lon,coords_lat
m.plot(x4,y4, linewidth = 2, color = 'k')


#UL bound lon+dx,lat+dy
x5,y5 = m(141.83,13.542)
m.scatter(x5,y5,s=200, color = 'k') #Far northern and western part of domain (tip of arrow w/ shear_dir 115 degrees)

# UR bound lon - dx, lat+dy #Far northern and eastern part of domain
#x6,y6 = m(140.92,15.29)
#m.scatter(x6,y6,s=200, color = 'm')



#UL bound lon+dx,lat+dy
#x5,y5 = m(140.083,14.87)
#m.scatter(x5,y5,s=200, color = 'k')

# UR bound lon - dx, lat+dy
#x6,y6 = m(140.92,14.87)
#m.scatter(x6,y6,s=200, color = 'm')


#Add the annuli
circle = Circle(xy = m(storm_lon,storm_lat), radius = 0.4054, fill = False, linewidth = 2)
plt.gca().add_patch(circle)

circle = Circle(xy = m(storm_lon,storm_lat), radius = 2.2522, fill = False, linewidth = 2)
plt.gca().add_patch(circle)


plt.title('Near-Surface KuPR 20230215-S011727-E024958', size = 28)
plt.show()
'''
