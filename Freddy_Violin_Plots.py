import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
import pandas as pd


files_dl  = glob.glob('0309_slope_warm_dl_eyewall*')
files_dr  = glob.glob('0309_slope_warm_dr_eyewall*')
files_ul  = glob.glob('0309_slope_warm_ul_eyewall*')
files_ur  = glob.glob('0309_slope_warm_ur_eyewall*')

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

print(np.nanmedian(dl_slope_warm))
print(np.nanmedian(dr_slope_warm))
print(np.nanmedian(ul_slope_warm))
print(np.nanmedian(ur_slope_warm))

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
plt.title(r'$\bf{i)}$ Slope of KuPR in Liquid Phase 03/09 1009 UTC (Svr TS)', size = title_font_size, y = 1.02)
ax.set_xticklabels(['DL','DR','UL','UR'],size = label_font_size)
ax.set_yticklabels(ax.get_yticks(), size = label_font_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'k', linewidth =  3.0, linestyle = '--')
plt.savefig('0309_warm.png', dpi = 300)
plt.show()
