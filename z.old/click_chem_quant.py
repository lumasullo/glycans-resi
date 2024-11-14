#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 16:38:31 2023

@author: masullo
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from scipy import signal
import os
import yaml

plt.close('all')

path = 'click chem quant/'
# path = 'click chem quant/DBCO/'

# filename = '231214_alkyne_preR1_R3-200pM-R2-100pM_1_MMStack_Pos0.ome_locs_render2_aligned_apicked.hdf5'
# title = 'Locs per pick histogram - R3'   
# binsize_exp = 1
# xmax = 100
# ymax = 70
# histcolor = '#2880C4'

filename = '231214_alkyne_postR1_R3-200pM-R2-100pM_1_MMStack_Pos0.ome_locs_render_aligned_apicked.hdf5'
title = 'Locs per pick histogram - R3'   
binsize_exp = 1
xmax = 100
ymax = 70
histcolor = '#2880C4'

# filename = '231213_azide__afterclick_overnight_R3-200pM-R2-100pM_1_MMStack_Pos0.ome_locs_render_aligned_apicked.hdf5'
# title = 'DBCO - Locs per pick histogram - R3'   
# binsize_exp = 2
# xmax = 100
# ymax = 70
# histcolor = '#2880C4'

# filename = '231209_01_azide-origami_20nm-grid_coloc-strand_after-click_R3-R5-200pM_1_MMStack_Pos0.ome_locs_render_aaligned_apicked_small.hdf5'
# title = 'Locs per pick histogram - R3'   
# binsize_exp = 2
# xmax = 100
# ymax = 70
# histcolor = '#2880C4'
 
# filename = '231209_01_azide-origami_20nm-grid_coloc-strand_R1-R5-200pM_1_MMStack_Pos0.ome_locs_render_aaligned_apicked.hdf5'
# title = 'Locs per pick histogram - R1'
# binsize_exp = 100
# xmax = 5000
# ymax = 80
# histcolor = '#97D8C4'

df = pd.read_hdf(path + filename, key='locs')

npicks = 928

threshold = 20
counter = 0
counter_zero = 0

nlocs_array = np.zeros(npicks)

for i in range(npicks):
    
    nlocs_array[i] = df[df['group'] == i].shape[0]
    
    if nlocs_array[i] > threshold:
    
        counter += 1
        
    elif nlocs_array[i] == 0:
        
        counter_zero += 1
        
    else:
        
        pass
    
    print(nlocs_array[i])
    
    
fig, ax = plt.subplots()
ax.set_title(title)

colors = ['#2880C4', '#97D8C4', '#F4B942', '#363636']

ax.set_xlim(0, xmax)
ax.set_ylim(0, ymax)

ax.set_xlabel('# of locks in pick')
ax.set_ylabel('Counts')

ax.set_box_aspect(1)

maxnum = nlocs_array.max()
bins_exp = np.arange(0, maxnum, binsize_exp)

counts, bin_edges, _ = ax.hist(nlocs_array, bins=bins_exp, 
                               edgecolor='black', linewidth=0.1, alpha=0.5, 
                               density=False, color=histcolor)

efficiency = counter / npicks

print('Click chemistry efficiency is ', np.around(efficiency, 3))
