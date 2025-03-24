#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 14:44:28 2022

@author: Luciano A. Masullo
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.neighbors import NearestNeighbors

plt.close('all')

# =============================================================================
# experimental parameters
# =============================================================================

# independent parameters

d = 2 # dimension of the simulation, d = 2 for 2D case, d = 3 for 3D
density_arr = np.linspace(1, 5000, 300) * 10**-6 # molecules per nm^2

σ_dnapaint_arr = np.array([3, 4]) # nm

width = 8e3 # width of the simulated area in nm
height = 8e3 # height of the simulated area in nm
distribution = 'uniform'

err_val = 0.10 # admitted frac of molecules closer than the resolution limit

# dependent parameters

resolution_arr = 4 * σ_dnapaint_arr # minimal distance between clusters to consider them resolvable
N = np.array(density_arr * width * height, dtype=int)

# =============================================================================
# simulate molecule positions
# =============================================================================

n_subres_frac_arr = np.zeros((len(density_arr), len(resolution_arr)))

for i, density in enumerate(density_arr):
    
    # calculate number of molecules
    N = int(density * width * height)

    # simulate positons for the molecules
    pos = np.zeros((N, d)) # initialize array of localizations
    
    if distribution == 'uniform':
    
        pos = np.array([np.random.uniform(0, width, N), 
                        np.random.uniform(0, height, N)]).T
    
    elif distribution == 'evenly spaced':
        
        wstep = width/np.sqrt(N)
        hstep = height/np.sqrt(N)
        pos = np.mgrid[0:width:wstep, 
                       0:height:hstep].reshape(2,-1).T
    
    
    nbrs = NearestNeighbors(n_neighbors=2).fit(pos) # find nearest neighbours
    _distances, _indices = nbrs.kneighbors(pos) # get distances and indices
    distances = _distances[:, 1] # get only the first neighbour distances
    
    for j, resolution in enumerate(resolution_arr):
    
        n_subres = len(distances[distances < resolution])
        n_subres_frac = n_subres/N
        
        n_subres_frac_arr[i, j] = n_subres_frac
        
# =============================================================================
# 1D plot frac of molecules unresolvable by RESI
# =============================================================================

fig2, ax2 = plt.subplots()

K = 6
density_arr_real = density_arr * K # RESI effectively decreases the density by a factor of K

for i in range(len(resolution_arr)):
    
    ax2.plot(density_arr_real * 10**6, n_subres_frac_arr[:, i], 
             label=str(np.around(resolution_arr[i], 0))+' nm')
    
ax2.set_xlabel('Density ($μm^{-2}$)')
ax2.set_ylabel('Fraction of subres molecules')
ax2.set_xlim(0, 5000)
ax2.tick_params(direction='in')

ax2.plot(density_arr * 10**6, np.ones(len(density_arr_real)) * err_val, 'k--')
    
ax2.legend()