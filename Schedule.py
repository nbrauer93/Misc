#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:40:23 2020

@author: noahbrauer
"""

import numpy as np
import matplotlib.pyplot as plt


months = np.arange(0,37,1)
x_ticks_labels = np.array(['General Exam', 'JPL Travel', 'Task 1', 'Task 2', 'Task 3', 'Dissertation Prep']).T[::-1]

##Create array with values to represent a color


color_array = np.zeros((6,36))*np.nan


color_array[5,0] = 10
color_array[4,1] = 30
color_array[3,0:7] = 50
color_array[2,6:24] = 35
color_array[1, 24:33] = 25
color_array[0,33:] = 70


fig = plt.figure(figsize = [15,10]) 
#plt.grid(color = 'k')


for i in range(len(months)):
    

    plt.axhline(i+0.49, color = 'k', linewidth = 0.75)
    plt.axvline(i+0.49, color = 'k', linewidth = 0.75)
    
    
im = plt.imshow(color_array, cmap = 'hsv')

plt.yticks(np.arange(0,7,1),x_ticks_labels, fontsize = 12)
plt.xticks(months, months, fontsize = 12)

plt.show()
