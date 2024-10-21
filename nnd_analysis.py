#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 12:26:20 2023

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

### Parameters ###

# parameters for the NND analysis
binsize = 1
maxdist = 1000

# display plots or not
plot_examples = False

fsize = (10, 10)
msize = 50
length = 1000 # nm, length of the display area for the graph

# pick info
# area_r = 40 # px
# px_size = 0.130 # μm
# pick_area = π * (area_r*px_size)**2
# pick_area = 1 * pick_area # quick fix for inhomogenous area

# pick_area = 75 * 75 # quick fix for full FOV

# pick_area = 18.62 # in μm^2

pick_area = 84.95 # in μm^2

# =============================================================================
# Load experimental data
# =============================================================================

# folder = '/Users/masullo/Documents/GitHub/glycans-resi/data/HMEC/ManNAz_june2024'
# file = 'mannaz_picked_resi.hdf5'
# file = 'mannaz_picked_resi_pick.hdf5'
# file = 'mannaz_picked_resi_pick_homogenous.hdf5'

# file = 'mannaz_picked_resi_pick_dense.hdf5'

# 240617_HMEC ManNAz
folder = r'/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/240617_HMEC/RESI_ManNAz/workflow_analysis/00_cluster_241009-1223/00_cluster_aggregation_241009-1223/04_save_datasets_aggregated'
file = r'target_picked'
title = '240617_ManNAz'

# folder = r'/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/240618_HMEC/RESI_GalNAz/workflow_analysis/00_cluster_241007-1052/00_cluster_aggregation_241007-1052/04_save_datasets_aggregated'
# file = r'target_picked'
# title = '240617_GalNAz'

# folder = 'data/HMEC'

# folder = 'data/MCF10A - MCF10AT'

# file = '231030_MCF10A-treated-TGFbeta_cell2_resi.hdf5'
# title = '231030_MCF10A-treated-TGFbeta'

# file = '231101_MCF10AT-treated-TGFbeta_cell1_resi.hdf5'
# title = '231101_MCF10AT-treated-TGFbeta'

# file = '231101_MCF10AT-treated-TGFbeta_cell2_resi.hdf5'
# title = '231101_MCF10AT-treated-TGFbeta'

# file = '231031_MCF10A-untreated_apicked_resi.hdf5'
# title = '231031_MCF10A-untreated'

# file = '231030_MCF10A-untreated_apicked_resi.hdf5'
# title = '231030_MCF10A-untreated'

# file = '231103_MCF10A-untreated_apicked_resi.hdf5'
# title = '231103_MCF10A-untreated'

# file = '231030_MCF10A-treated-TGFbeta_apicked_resi.hdf5'
# title = '231030_MCF10A-treated-TGFbeta'

# file = '231031_MCF10A-treated-TGFbeta_apicked_resi.hdf5'
# title = '231031_MCF10A-treated-TGFbeta'

# file = '231101_MCF10AT-untreated_apicked_resi.hdf5'
# title = '231101_MCF10AT-untreated'

# file = '231102_MCF10AT-untreated_apicked_resi.hdf5'
# title = '231102_MCF10AT-untreated'

# file = '231103_MCF10AT-treated-TGFbeta_apicked_resi.hdf5'
# title = '231103_MCF10AT-treated-TGFbeta'

# file = '231106_MCF10AT-treated-TGFbeta_apicked_resi.hdf5'
# title = '231106_MCF10AT-treated-TGFbeta'

# file = '231101_MCF10AT-treated-TGFbeta_apicked_resi.hdf5'
# title = '231101_MCF10AT-treated-TGFbeta'

# file = 'ManNAz_cell1_resi.hdf5'
# title = 'HMEC_ManNAz'

# file = '230619_ManNAz_RESI.hdf5'
# title = '230619_HMEC_ManNAz - full FOV'

# file = '230619_ManNAz_RESI_picked.hdf5'
# title = '230619_HMEC_ManNAz'

# file = '230620_ManNAz_RESI_picked.hdf5'
# title = '230620_HMEC_ManNAz'

# file = '230619_ManNAz_RESI_picked_highdensity.hdf5'
# title = '230619_HMEC_ManNAz - highdensity area'

filename_0 = folder + '/' + file + '.hdf5'

df_0 = pd.read_hdf(filename_0, key = 'locs')

x_0 = df_0.x*130 # in nm
y_0 = df_0.y*130 # in nm

pos_exp_0 = np.array([x_0, y_0]).T

n_sugars = pos_exp_0.shape[0]
obs_density = n_sugars / pick_area # molec / μm^2

print(title)
print('N sugars = ', n_sugars)
print('Density = ' + str(int(np.around(obs_density,0))) + '/μm^2')

# =============================================================================
# NN calculation for experimental data
# =============================================================================

# find nearest neighbours from biomolecule 0 to biomolecule 1 (exp data), in this case is the same biomolecule (e.g. ManNAz, GalNAz)
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

ax0.set_title(title)

colors = ['#2880C4', '#97D8C4', '#F4B942', '#363636']

ax0.set_xlim(0, 100)
ax0.set_ylim(0, 0.06)

ax0.set_xlabel('K-th NND (nm)')
ax0.set_ylabel('Counts')

ax0.set_box_aspect(1)

binsize_exp = 0.5
bins_exp = np.arange(0, maxdist, binsize_exp)

for i in range(4):
    
    
    key = str(i+1) + 'nn'

    counts, bin_edges, _ = ax0.hist(distances_exp[key], bins=bins_exp, 
                                    edgecolor='black', linewidth=0.1, alpha=0.5, 
                                    density=True, color=colors[i], label=key)

# ax0.legend()

# =============================================================================
# Simulation of Complete Spatial Randomness (CSR)
# =============================================================================

N = int(1e6)
width = np.sqrt(N/obs_density) 
height = width

# simple CSR model
pos_sim = np.array([np.random.uniform(0, width, N), 
                    np.random.uniform(0, height, N)]).T

pos_sim = pos_sim * 1000 # in nm

# =============================================================================
# NN calculation for simulated data
# =============================================================================

# find nearest neighbours from protein 0 to protein 1 (exp data)
nbrs = NearestNeighbors(n_neighbors=5).fit(pos_sim) 
_distances_sim, _indices_sim = nbrs.kneighbors(pos_sim) # get distances and indices

distances_sim = {}
freqs_sim = {}
bin_centers_sim = {}

binsize_sim = 1
bins_sim = np.arange(0, maxdist, binsize_sim)

# calculate K-th NND
for i in range(4):
    
    key = str(i+1) + 'nn'
    distances_sim[key] = _distances_sim[:, i+1] # get the first neighbour distances (0 for hetero, 1 for homo in the second coord)

    bins_sim = np.arange(0, maxdist, binsize_sim)
    freqs_sim[key], bin_edges_sim = np.histogram(distances_sim[key], bins=bins, density=True)

    bin_centers_sim[key] = (bin_edges_sim[:-1] + bin_edges_sim[1:])/2
    
# plot simulation    
for i in range(4):
    
    key = str(i+1) + 'nn'
    ax0.plot(bin_centers_sim[key], freqs_sim[key], color=colors[i], linewidth=2)

# ax0.set_xlim(0, 200)
# ax0.set_ylim(0, freqs_sim['1nn'].max()*1.1)

ax0.set_xlabel('K-th NND (nm)')
ax0.set_ylabel('Frequency')

ax0.set_box_aspect(1)



