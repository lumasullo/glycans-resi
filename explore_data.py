#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:31:09 2023

@author: Luciano A. Masullo
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from datetime import datetime
import configparser
import time 

π = np.pi
plt.close('all')

### Misc parameters ###

# parameters for the NND analysis
binsize = 2
maxdist = 1000

# display plots or not
plot_examples = False

fsize = (10, 10)
msize = 50
length = 1000 # nm, length of the display area for the graph

# =============================================================================
# Load experimental data
# =============================================================================

folder = 'data/230619'
file = 'ManNAz_RESI.hdf5'
# file = 'GalNAz_RESI.hdf5'
# file = 'ManNAz_R1_250pM_pcapcdtx_561_35mW_1_MMStack_Pos0.ome_locs_render_RCC500_undriftfrompicks_without_AuNPs_cluster_centers.hdf5'
# file = 'ManNAz_R2_150pM_pcapcdtx_561_35mW_1_MMStack_Pos0.ome_locs_render_RCC500_undriftfrompicks_aligned_without_AuNPs_cluster_centers.hdf5'
# file = 'ManNAz_R3_250pM_pcapcdtx_561_35mW_1_MMStack_Pos0.ome_locs_render_RCC500_undriftbypicks_aligned_without_AuNPs_cluster_centers.hdf5'
# file = 'ManNAz_R4_250pM_pcapcdtx_561_35mW_1_MMStack_Pos0.ome_locs_render_RCC500_undriftfrompicks_aligned_without_AuNPs_cluster_centers.hdf5'
# file = 'ManNAz_R5_500pM_pcapcdtx_561_35mW_1_MMStack_Pos0.ome_locs_render_RCC500_undriftbypicks_aligned_without_AuNPs_cluster_centers.hdf5'
# file = 'ManNAz_R6_250pM_pcapcdtx_561_35mW_1_MMStack_Pos0.ome_locs_render_RCC500_undriftfrompicks_aligned_without_AuNPs_cluster_centers.hdf5'

filename_0 = folder + '/' + file 

df_0 = pd.read_hdf(filename_0, key = 'locs')

x_0 = df_0.x*130
y_0 = df_0.y*130

pos_exp_0 = np.array([x_0, y_0]).T

# =============================================================================
# NN calculation for experimental data
# =============================================================================

# find nearest neighbours from protein 0 to protein 1 (exp data)
nbrs = NearestNeighbors(n_neighbors=5).fit(pos_exp_0) 
_distances_exp, _indices_exp = nbrs.kneighbors(pos_exp_0) # get distances and indices

distances_exp = {}
freqs_exp = {}
bin_centers_exp = {}

for i in range(4):
    
    key = str(i+1) + 'nn'
    distances_exp[key] = _distances_exp[:, i+1] # get the first neighbour distances (0 for hetero, 1 for homo in the second coord)

    bins = np.arange(0, maxdist, binsize)
    freqs_exp[key], bin_edges = np.histogram(distances_exp[key], bins=bins, density=True)

    bin_centers_exp[key] = (bin_edges[:-1] + bin_edges[1:])/2
    
# =============================================================================
# Plot exp data
# =============================================================================

fig0, ax0 = plt.subplots()

colors = ['#2880C4', '#97D8C4', '#F4B942', '#363636']

ax0.set_xlim(0, 200)
# ax0.set_ylim(0, freqs_sim['1nn'].max()*1.1)

ax0.set_xlabel('K-th NND (nm)')
ax0.set_ylabel('Counts')

ax0.set_box_aspect(1)

binsize = 2
bins = np.arange(0, maxdist, binsize)

for i in range(4):
    
    
    key = str(i+1) + 'nn'

    counts, bin_edges, _ = ax0.hist(distances_exp[key], bins=bins, 
                                    edgecolor='black', linewidth=0.1, alpha=0.5, 
                                    density=True, color=colors[i], label=key)

ax0.legend()
