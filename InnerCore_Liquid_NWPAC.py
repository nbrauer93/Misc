#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 17:05:56 2022

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
import pandas as pd

#Import file list for each shear-relative quadrant

files_dl  = glob.glob('*_slope_warm_dl_inner*')
files_dr  = glob.glob('*_slope_warm_dr_inner*')
files_ul  = glob.glob('*_slope_warm_ul_inner*')
files_ur  = glob.glob('*_slope_warm_ur_inner*')

#Now read in each file


dl_slope = []
dr_slope = []
ul_slope = []
ur_slope = []

for i in range(len(files_dl)):
    
    dl_open = np.load(files_dl[i], allow_pickle = False)
    dl_slope.append(dl_open)
    
    dr_open = np.load(files_dr[i], allow_pickle = False)
    dr_slope.append(dr_open)
    
    ul_open = np.load(files_ul[i], allow_pickle = False)
    ul_slope.append(ul_open)
    
    ur_open = np.load(files_ur[i], allow_pickle = False)
    ur_slope.append(ur_open)


#Flatten lists to arrays    
    
dl_slope_warm = np.concatenate(dl_slope).ravel()    
dr_slope_warm = np.concatenate(dr_slope).ravel()  
ul_slope_warm = np.concatenate(ul_slope).ravel()  
ur_slope_warm = np.concatenate(ur_slope).ravel()  





#Now plot as box and whisker plots



data_array = [dl_slope_warm,dr_slope_warm,ul_slope_warm,ur_slope_warm]


label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array)
#
plt.ylim(-2,2)
plt.ylabel('[dBZ/km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Slope of KuPR in Outer Annulus in Liquid Phase (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
#ax.set_yticklabels(np.arange(-2.25,2.25,step = 0.25),size = label_font_size)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)

plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'k', linewidth =  3.0, linestyle = '--')

#%%

iles_dl  = glob.glob('*_slope_warm_dl_eyewall*')
files_dr  = glob.glob('*_slope_warm_dr_eyewall*')
files_ul  = glob.glob('*_slope_warm_ul_eyewall*')
files_ur  = glob.glob('*_slope_warm_ur_eyewall*')

#Now read in each file


dl_slope = []
dr_slope = []
ul_slope = []
ur_slope = []

for i in range(len(files_dl)):

    dl_open = np.load(files_dl[i], allow_pickle = False)
    dl_slope.append(dl_open)

    dr_open = np.load(files_dr[i], allow_pickle = False)
    dr_slope.append(dr_open)

    ul_open = np.load(files_ul[i], allow_pickle = False)
    ul_slope.append(ul_open)

    ur_open = np.load(files_ur[i], allow_pickle = False)
    ur_slope.append(ur_open)


#Flatten lists to arrays

dl_slope_warm = np.concatenate(dl_slope).ravel()
dr_slope_warm = np.concatenate(dr_slope).ravel()
ul_slope_warm = np.concatenate(ul_slope).ravel()
ur_slope_warm = np.concatenate(ur_slope).ravel()





#Now plot as box and whisker plots



data_array = [dl_slope_warm,dr_slope_warm,ul_slope_warm,ur_slope_warm]


label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.ylim(-2,2)
plt.ylabel('[dBZ/km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Slope of KuPR in Inner Annulus in Liquid Phase (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'k', linewidth =  3.0, linestyle = '--')




#%%


#Now do the exact same thing for the ice phase:
    
files_dl_ice  = glob.glob('*_slope_ice_dl_inner*')
files_dr_ice  = glob.glob('*_slope_ice_dr_inner*')
files_ul_ice  = glob.glob('*_slope_ice_ul_inner*')
files_ur_ice  = glob.glob('*_slope_ice_ur_inner*')

#Now read in each file


dl_slope_ice = []
dr_slope_ice = []
ul_slope_ice = []
ur_slope_ice = []

for i in range(len(files_dl_ice)):
    
    dl_open = np.load(files_dl_ice[i], allow_pickle = False)
    dl_slope_ice.append(dl_open)
    
    dr_open = np.load(files_dr_ice[i], allow_pickle = False)
    dr_slope_ice.append(dr_open)
    
    ul_open = np.load(files_ul_ice[i], allow_pickle = False)
    ul_slope_ice.append(ul_open)
    
    ur_open = np.load(files_ur_ice[i], allow_pickle = False)
    ur_slope_ice.append(ur_open)


#Flatten lists to arrays    
    
dl_slope_ice = np.concatenate(dl_slope_ice).ravel()    
dr_slope_ice = np.concatenate(dr_slope_ice).ravel()  
ul_slope_ice = np.concatenate(ul_slope_ice).ravel()  
ur_slope_ice = np.concatenate(ur_slope_ice).ravel()  





#Now plot as box and whisker plots



data_array_ice = [dl_slope_ice,dr_slope_ice,ul_slope_ice,ur_slope_ice]


label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array_ice)
#

plt.ylabel('[dBZ/km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Slope of KuPR in Outer Annulus in Ice Phase (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
#ax.set_yticklabels(np.arange(-2.25,2.25,step = 0.25),size = label_font_size)
plt.ylim(-2,2) 
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)

plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'k', linewidth =  3.0, linestyle = '--')

#%%


files_dl_ice  = glob.glob('*_slope_ice_dl_eyewall*')
files_dr_ice  = glob.glob('*_slope_ice_dr_eyewall*')
files_ul_ice  = glob.glob('*_slope_ice_ul_eyewall*')
files_ur_ice  = glob.glob('*_slope_ice_ur_eyewall*')

#Now read in each file


dl_slope_ice = []
dr_slope_ice = []
ul_slope_ice = []
ur_slope_ice = []

for i in range(len(files_dl_ice)):
    
    dl_open = np.load(files_dl_ice[i], allow_pickle = False)
    dl_slope_ice.append(dl_open)
    
    dr_open = np.load(files_dr_ice[i], allow_pickle = False)
    dr_slope_ice.append(dr_open)
    
    ul_open = np.load(files_ul_ice[i], allow_pickle = False)
    ul_slope_ice.append(ul_open)
    
    ur_open = np.load(files_ur_ice[i], allow_pickle = False)
    ur_slope_ice.append(ur_open)


#Flatten lists to arrays    
    
dl_slope_ice = np.concatenate(dl_slope_ice).ravel()    
dr_slope_ice = np.concatenate(dr_slope_ice).ravel()  
ul_slope_ice = np.concatenate(ul_slope_ice).ravel()  
ur_slope_ice = np.concatenate(ur_slope_ice).ravel()  





#Now plot as box and whisker plots



data_array_ice = [dl_slope_ice,dr_slope_ice,ul_slope_ice,ur_slope_ice]


label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array_ice)
#

plt.ylabel('[dBZ/km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Slope of KuPR in Inner Annulus in Ice Phase (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
#ax.set_yticklabels(np.arange(-2.25,2.25,step = 0.25),size = label_font_size)
plt.ylim(-2,2) 
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)

plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'k', linewidth =  3.0, linestyle = '--')
   

#%%

#Now do echo tops

files_dl_echo  = glob.glob('*_top_dl_inner*')
files_dr_echo  = glob.glob('*_top_dr_inner*')
files_ul_echo  = glob.glob('*_top_ul_inner*')
files_ur_echo  = glob.glob('*_top_ur_inner*')

#Now read in each file


dl_top = []
dr_top = []
ul_top = []
ur_top = []

for i in range(len(files_dl_echo)):
    
    dl_open = np.load(files_dl_echo[i], allow_pickle = False)
    dl_top.append(dl_open)
    
    dr_open = np.load(files_dr_echo[i], allow_pickle = False)
    dr_top.append(dr_open)
    
    ul_open = np.load(files_ul_echo[i], allow_pickle = False)
    ul_top.append(ul_open)
    
    ur_open = np.load(files_ur_echo[i], allow_pickle = False)
    ur_top.append(ur_open)


#Flatten lists to arrays    
    
dl_top = np.concatenate(dl_top).ravel()    
dr_top = np.concatenate(dr_top).ravel()  
ul_top = np.concatenate(ul_top).ravel()  
ur_top = np.concatenate(ur_top).ravel()  





#Now plot as box and whisker plots



data_array_top = [dl_top,dr_top,ul_top,ur_top]


label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array_top)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.ylim(0,20)
plt.ylabel('[km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Echo Top Height in Outer Annulus (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
    

#%%



#Now do echo tops

files_dl_echo  = glob.glob('*_top_dl_eyewall*')
files_dr_echo  = glob.glob('*_top_dr_eyewall*')
files_ul_echo  = glob.glob('*_top_ul_eyewall*')
files_ur_echo  = glob.glob('*_top_ur_eyewall*')

#Now read in each file


dl_top = []
dr_top = []
ul_top = []
ur_top = []

for i in range(len(files_dl_echo)):
    
    dl_open = np.load(files_dl_echo[i], allow_pickle = False)
    dl_top.append(dl_open)
    
    dr_open = np.load(files_dr_echo[i], allow_pickle = False)
    dr_top.append(dr_open)
    
    ul_open = np.load(files_ul_echo[i], allow_pickle = False)
    ul_top.append(ul_open)
    
    ur_open = np.load(files_ur_echo[i], allow_pickle = False)
    ur_top.append(ur_open)


#Flatten lists to arrays    
    
dl_top = np.concatenate(dl_top).ravel()    
dr_top = np.concatenate(dr_top).ravel()  
ul_top = np.concatenate(ul_top).ravel()  
ur_top = np.concatenate(ur_top).ravel()  





#Now plot as box and whisker plots



data_array_top = [dl_top,dr_top,ul_top,ur_top]


label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array_top)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.ylim(0,20)
plt.ylabel('[km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Echo Top Height in Inner Annulus (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
    



#%%


#Plot uncorrected KuPR profiles



    
files_dl_ice  = glob.glob('*_slope_uncor_ice_dl_eyewall*')
files_dr_ice  = glob.glob('*_slope_uncor_ice_dr_eyewall*')
files_ul_ice  = glob.glob('*_slope_uncor_ice_ul_eyewall*')
files_ur_ice  = glob.glob('*_slope_uncor_ice_ur_eyewall*')

#Now read in each file


dl_slope_ice = []
dr_slope_ice = []
ul_slope_ice = []
ur_slope_ice = []

for i in range(len(files_dl_ice)):
    
    dl_open = np.load(files_dl_ice[i], allow_pickle = False)
    dl_slope_ice.append(dl_open)
    
    dr_open = np.load(files_dr_ice[i], allow_pickle = False)
    dr_slope_ice.append(dr_open)
    
    ul_open = np.load(files_ul_ice[i], allow_pickle = False)
    ul_slope_ice.append(ul_open)
    
    ur_open = np.load(files_ur_ice[i], allow_pickle = False)
    ur_slope_ice.append(ur_open)


#Flatten lists to arrays    
    
dl_slope_ice = np.concatenate(dl_slope_ice).ravel()    
dr_slope_ice = np.concatenate(dr_slope_ice).ravel()  
ul_slope_ice = np.concatenate(ul_slope_ice).ravel()  
ur_slope_ice = np.concatenate(ur_slope_ice).ravel()  


 
data_array_ice = [dl_slope_ice,dr_slope_ice,ul_slope_ice,ur_slope_ice]

label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array_ice)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.ylim(-2,2)
plt.ylabel('[dBZ/km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Slope of KuPR (Uncorr) in Eyewall in Ice Phase (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'k', linewidth =  3.0, linestyle = '--')

#%%



    
files_dl_ice  = glob.glob('*_slope_uncor_ice_dl_inner*')
files_dr_ice  = glob.glob('*_slope_uncor_ice_dr_inner*')
files_ul_ice  = glob.glob('*_slope_uncor_ice_ul_inner*')
files_ur_ice  = glob.glob('*_slope_uncor_ice_ur_inner*')

#Now read in each file


dl_slope_ice = []
dr_slope_ice = []
ul_slope_ice = []
ur_slope_ice = []

for i in range(len(files_dl_ice)):
    
    dl_open = np.load(files_dl_ice[i], allow_pickle = False)
    dl_slope_ice.append(dl_open)
    
    dr_open = np.load(files_dr_ice[i], allow_pickle = False)
    dr_slope_ice.append(dr_open)
    
    ul_open = np.load(files_ul_ice[i], allow_pickle = False)
    ul_slope_ice.append(ul_open)
    
    ur_open = np.load(files_ur_ice[i], allow_pickle = False)
    ur_slope_ice.append(ur_open)


#Flatten lists to arrays    
    
dl_slope_ice = np.concatenate(dl_slope_ice).ravel()    
dr_slope_ice = np.concatenate(dr_slope_ice).ravel()  
ul_slope_ice = np.concatenate(ul_slope_ice).ravel()  
ur_slope_ice = np.concatenate(ur_slope_ice).ravel()  


 
data_array_ice = [dl_slope_ice,dr_slope_ice,ul_slope_ice,ur_slope_ice]

label_font_size = 20
title_font_size = 20

plt.figure(figsize=(10,10))
ax = sns.violinplot(data=data_array_ice)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.ylim(-2,2)
plt.ylabel('[dBZ/km]', fontsize = label_font_size)
plt.xlabel('850-200 hPa Shear-Relative Quadrant', fontsize = label_font_size)
plt.title('Slope of KuPR (Uncorr) in Inner Core in Ice Phase (NWPAC)', size = title_font_size)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'k', linewidth =  3.0, linestyle = '--')