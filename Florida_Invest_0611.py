import h5py
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj
from mpl_toolkits.basemap import Basemap


file = '2A.GPM.DPR.V9-20240130.20240611-S144615-E161927.058423.V07C.HDF5'


#Read in file, extract lat, lons, PSD

DPR = h5py.File(file, 'r')

#nscan is the nummber of scans in each granule
#nray is the number of angle bins in each  scan; think of these as footprint scans (5 km in diameter for each footprint)
#nbin is the number of range bins in each ray (angle bins)
#nDSD: Parameters are log_{10}(N_{w}) (number concentration) and D0 (mean drop diameter)


lat = DPR['FS']['Latitude'][:,:]
lon = DPR['FS']['Longitude'][:,:]
z = DPR['FS']['SLV']['zFactorFinalESurface'][:,:,0] #nscan x nray (7934,49)
dsd = DPR['FS']['SLV']['paramDSD'][:]
freezing = DPR['FS']['VER']['heightZeroDeg'][:]

#%%

###Define cross-section length

ind1 = np.where((lon[:,0]>=-78.5)) #Where are the DPR lons >= -100
ind2 = np.where((lon[:,0])<=-77.5) #Where are lons <= -85
ind3 = np.intersect1d(ind1,ind2) #Both conditions need to be satisfied here


#Change the ray to change cross-section location (i.e. the "27" in this case)


###Setup to 2D grid for plotting

x = 2.* 17 #48 degrees (from -17 to 17)
re = 6378. #radius of the earth
theta = -1*(x/2.) + (x/48.)*np.arange(0,49) #Split into equal degrees (from -17 to 17)
theta2  = np.ones(theta.shape[0]+1)*np.nan #Define an empty array (NaNs) with same shape as ray dimension
theta = theta - 0.70833333/2. #Shift degrees for plotting pruposes
theta2[:-1] = theta #remove last array dimension in theta so python won't get angry about shape
theta2[-1] = theta[-1] + 0.70833333
theta = theta2*(np.pi/180.) #Convert from degrees to radians

prh = np.ones((177,50))*np.nan #Define empty grid

for i in range(prh.shape[0]): #Loop over each range gate
    for j in range(prh.shape[1]): #Loop over each scan
            a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #Orbit height of 407 km

            prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a)

h2 = prh #where prh is the (range bins,ray)-space
h3 =np.ones((177,50))*np.nan

for i in range(h3.shape[1]):
    h3[:,i] = h2[::-1,i] #This reverses the vertical dimension so as indices increase, height increases


ku = DPR['FS']['SLV']['zFactorFinal'][ind3,:,:] #Read in ku-band reflectivity; nscan x nray (554,49,176)
n0 = dsd[ind3,:,:,0]/10 #Read in the number concentration
d0 = dsd[ind3,:,:,1] #Read in the mean drop diameter  #Both have dimensions nscan x nray x nbin (554,49,176)
zeroDeg = freezing[ind3,:]

#Cut all parameters so they are at same ray as above
#45 is default

ray = 31

ku = ku[:,ray,:]
n0 = n0[:,ray,:]
d0 = d0[:,ray,:]
zero_deg_isotherm = zeroDeg[:,ray]/1000 #Convert from meters to km

#Take lats and lons along same ray
lons = DPR['FS']['Longitude'][ind3,ray]
lats = DPR['FS']['Latitude'][ind3,ray]




#Choose a starting point, then calculate distance
lat0 = lats[0]
lon0 = lons[0]


p = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0) #Define a projection and set starting lat an lon to same point as above


#%%

#Define a 2D array for plotting purposes

lat_3d = np.ones(ku.shape)*np.nan
lon_3d = np.ones(ku.shape)*np.nan

for i in range(ku.shape[0]):
    lat_3d[i,:] = lats[i]
    lon_3d[i,:] = lons[i]


x,y = p(lon_3d,lat_3d) #Now convert degrees to distance (in km)
R_gpm = np.sqrt(x**2 + y**2)*np.sign(x) #Keeps sign of number; converts to radial distance

#Reverse range gate order for all parameters

ku = ku[:,::-1,0]
n0 = n0[:,::-1]
d0 = d0[:,::-1]



ku = np.ma.masked_where(ku<=12, ku) #Mask all the bad points in ku data
y = np.ones([ku.shape[0], ku.shape[1]]) #Define an empty array

#Define the appropriate range bins
h4 = h3[:,ray] #This number MUST match the same ray being used
for i in range(y.shape[1]):
    y[:,i] = h4[i]

#%%

print(len(lons))


#Determine median value of ind3 (intersection between longitudes)
cross_track_index = int(len(ind3)/2)+ind3[0]



#Let's plot a map with the center point location for our cross-sections

#Mask near surface reflectivity values

z = np.ma.masked_where(z<=12, z)

len_lons = len(lons)
len_lats = len(lats)

cmin = 12.; cmax = 70.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)

plt.figure(figsize=(10,10))

#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-82,-76]); ylim = np.array([24,30])
#xlim = np.array([35,55]); ylim = np.array([-30,-10])



m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawstates(), m.drawcountries()
cs = m.contourf(lon.T,lat.T,z.T,clevs,cmap='turbo',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lons[0:len_lons],lats[0:len_lats], 'k--', linewidth = 4)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()
parallels = np.arange(-60,60,step = 2)

m.drawparallels(parallels, labels = [True, False, False, False], fontsize = 20)


meridians = np.arange(0,360, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True],fontsize = 20)
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[dBZ]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    #cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)

plt.title('06/12 1549 UTC KuPR (Attenuation-Corrected)', size = 20, y = 1.04)

#plt.show()



#Remove the values less than or equal to zero

n0[n0<=0] = np.nan
d0[d0<=0] = np.nan



#%%

#Now we plot! #N-S (along track first)


#Plot mean drop size first

dist = np.array([0,300])
'''
x_loc = 500
y_loc = 14.5
label = '@NOAABrauer'
'''
plt.figure(figsize=(10,10))

vmax = 3
vmin = 0

R_min = R_gpm.min()
R_max = R_gpm.max()

R_gpm = R_gpm[:,:,0]

label_size = 20

cmin =0.; cmax = 3.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, d0, cmap='turbo',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'06/12 1549 UTC $D_{M}$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(0,14)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[mm]',size = label_size)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.clim(0,3)


#plt.text(x_loc, y_loc, label)


plt.show()



####Number concentration (liquid water content)

plt.figure(figsize=(10,10))

vmax = 6
vmin = 1

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =1.; cmax = 6.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, n0, cmap='turbo',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'06/12 1549 UTC $log_{10 }(N_{W})$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(0,14)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[$m^{-3}$ $mm^{-1}$]',size = label_size)
plt.clim(1,6)
#plt.text(x_loc, y_loc, label)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.show()


###And lastly Ku-band


plt.figure(figsize=(10,10))

vmax = 60
vmin =12

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, ku, cmap='turbo',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'06/12 1549 UTC KuPR ', size = 20)
#plt.xlim(300,450)
plt.xlim(dist[0], dist[1])
plt.ylim(0,14)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.clim(12,60)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

#plt.text(x_loc, y_loc, label)
plt.show()


#%%
#Now plot a map of slope:


ku_slope = DPR['FS']['SLV']['zFactorFinal'][:,:,:,0]
freezing_slope =  DPR['FS']['VER']['heightZeroDeg'][:]
lat_slope = DPR['FS']['Latitude'][:]
lon_slope = DPR['FS']['Longitude'][:]


#Reshape DPR lat/lon arrays to 1D space:



I,J,K = ku_slope.shape
ku_reshape = ku_slope.reshape(I*J,K, order = 'F')
freezing_reshape = freezing_slope.reshape(I*J, order = 'F')/1000
latitude = lat_slope.reshape(I*J, order = 'F')
longitude = lon_slope.reshape(I*J, order = 'F')


ku_reshape = ku_reshape[:,::-1]

x_bins = np.arange(10,60,1)
y_bins = np.arange(0,21.988032,0.124932)


#Normalize height w.r.t 0degC isotherm

norm_alt =[]

for i in range(len(freezing_reshape)):

    norm_altitude = y_bins - freezing_reshape[i]
    norm_alt.append(norm_altitude)

#Now do all the masking in order to perform linear regression

norm_alt = np.asarray(norm_alt)


norm_altitude = norm_alt.copy()
#norm_altitude[(norm_altitude>-1)|(norm_altitude<-3)] = -9999
#norm_altitude = np.ma.masked_where(norm_altitude<-10,norm_altitude)
#norm_altitude = np.ma.masked_where(norm_altitude<-3,norm_altitude)


ku_reshaped = ku_reshape.copy()
#ku_reshaped[(norm_altitude>-1)|(norm_altitude<-3)] = -9999
ku_reshaped[(norm_altitude>-1)|(norm_altitude<-3)] = np.nan
#ku_reshaped = np.ma.masked_where(ku_reshaped<0,ku_reshaped)



'''
'''
#Now compute the linear regression
#The bug is in here:
from scipy import stats
#from sklearn.linear_model import LinearRegression

slope_warm = []

for i in range(ku_reshaped.shape[0]):



    if len(ku_reshaped[i][~np.isnan(ku_reshaped[i])]) == 0:
        slope_warm.append(np.nan)
        continue

    #Slice all points of ku_reshaped to remove any values l

    sliced_norm_altitude = norm_altitude[i][(norm_altitude[i] <= -1) & (norm_altitude[i] >= -3)]
    sliced_ku_reshaped = ku_reshaped[i][(norm_altitude[i] <= -1) & (norm_altitude[i] >= -3)]

    filtered_norm_altitude = sliced_norm_altitude[~np.isnan(sliced_ku_reshaped)]
    filtered_ku_reshaped = sliced_ku_reshaped[~np.isnan(sliced_ku_reshaped)]

    if len(filtered_ku_reshaped) <= 5:
        slope_warm.append(np.nan)
        continue


    slope = stats.linregress(filtered_ku_reshaped,filtered_norm_altitude)[0]
    #print(slope)
    #model = LinearRegression().fit(ku_reshaped[i,:],norm_altitude[i,:])
    #slope = model.score(ku_reshaped[i,:],norm_altitude[i,:])
    slope_warm.append(slope)





slope_warm = np.asarray(slope_warm)

slope_liquid = slope_warm.reshape(I,J, order = 'F')

slope_liquid[(slope_liquid<-2)| (slope_liquid>2)] = np.nan



#%%
#Let's plot on a grid now

from mpl_toolkits.basemap import Basemap

cmin = -1.25; cmax = 1.25; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)

plt.figure(figsize=(10,10))





label_size = 22

xlim = np.array([-97,-91]); ylim = np.array([26,32])
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,slope_liquid,clevs,cmap='bwr',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
parallels = np.arange(-90,90,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False], fontsize = 20)

meridians = np.arange(0,360, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True], fontsize = 20)

cbar = plt.colorbar(fraction=0.035, pad=0.04)
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dB/km]',size = label_size)
plt.title('Vertical Slope of KuPR in Liquid Phase 06/12 1549 UTC', size = 20, y = 1.04)
plt.show()
