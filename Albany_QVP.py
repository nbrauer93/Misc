import matplotlib.pyplot as plt
import numpy as np
import pyart
import cartopy.crs as ccrs
import glob
#from tqdm import tqdm
from scipy import stats

radar_file = 'KENX20230314_000156_V06'



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



#Now do this for multiple files



files = glob.glob('KENX*')

#Loop through file list to get times; Parse out of file radar_name
#Sort file names by time first:

order_index = np.argsort(files)
files_ordered = [files[i] for i in order_index]

times = []
qvps_z = np.ones((237,1832))*np.nan
qvps_zdr = np.ones((237,1832))*np.nan
qvps_rhohv = np.ones((237,1832))*np.nan
#qvps_kdp = np.ones((107,1832))*np.nan

for i in range(len(files_ordered)):

    qvps = quasi_vertical_profile(files_ordered[i])

    qvps_z[i,:] = qvps['reflectivity']
    qvps_zdr[i,:] = qvps['differential_reflectivity']
    qvps_rhohv[i,:] = qvps['cross_correlation_ratio']
    #qvps_kdp[i,:] = pyart.retrieve.kdp_maesaka(qvps['differential_phase'])[0]

    parsed_time = files_ordered[i].split("_")[1]
    print(parsed_time)
    times.append(parsed_time)

#print(qvps_z.shape) # 107 x 1832
#print(np.nanmax(qvps_z))  # 44.495 dBZ

times = np.asarray(times)



#%%

qvps_z[qvps_z<5] = np.nan
qvps_zdr[qvps_zdr<0] = np.nan
qvps_rhohv[qvps_rhohv<0.8] = np.nan
qvps_rhohv[qvps_rhohv>1] = np.nan


#Now plot:

font_size = 16
title_size = 20

plt.figure(figsize=(10,10))

qvps_z[qvps_z<0] = np.nan
qvps_zdr[qvps_zdr<0] = np.nan

plt.pcolormesh(times,qvp['height']/1000.0,qvps_z.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size= font_size)
plt.ylabel('Height (km)', size = font_size)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = '[dBZ]',size = font_size)
#plt.colorbar()
plt.ylim(0,8)

plt.title(r'KENX Quasi-Vertical Profile 3/14 $Z_{H}$', size = title_size)
x = [0,59,118,177,236]
labels = np.array(['0000','0600','1200','1800','0000'])
plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()


plt.figure(figsize=(10,10))

plt.pcolormesh(times,qvp['height']/1000.0,qvps_zdr.T, cmap = 'pyart_RefDiff')
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('Height (km)', size = font_size)
plt.ylim(0,8)
plt.clim(0,3)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = '[dB]',size = font_size)
plt.title(r'KENX Quasi-Vertical Profile 3/14 $Z_{DR}$', size = title_size)

plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()

plt.figure(figsize=(10,10))

plt.pcolormesh(times,qvp['height']/1000.0,qvps_rhohv.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('Height (km)', size = font_size)
plt.ylim(0,8)
plt.clim(0.9,1.)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = r'[$\rho_{HV}$]',size = font_size)
plt.title(r'KENX Quasi-Vertical Profile 3/14 $\rho_{HV}$', size = title_size)

plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()
