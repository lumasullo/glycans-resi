#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:41:10 2024

@author: masullo
"""

import numpy as np
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import configparser
from timeit import default_timer as timer
from scipy.stats import gaussian_kde

from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import pdist


abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

from Functions import io
from Functions import dbscan

plt.close('all')

"""
===============================================================================
Load experimental data
===============================================================================
"""

# # 240617_HMEC ManNAz
# path = r'/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/240617_HMEC/RESI_ManNAz/workflow_analysis/00_cluster_241009-1223/00_cluster_aggregation_241009-1223/04_save_datasets_aggregated/'
# filename = r'target_picked.hdf5'

# #240618_HMEC GalNAz
path = r'/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/240618_HMEC/RESI_GalNAz/workflow_analysis/00_cluster_241007-1052/00_cluster_aggregation_241007-1052/04_save_datasets_aggregated/'
filename = r'target_picked.hdf5'

# path = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/Glycosylated spherical domains HMECs/Spherical clusters/240617HMECmannaz/1/'
# path = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/Homogenous areas/GalNAz/1/'

# path = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/ManNAz/240617_ManNAz/1/'
# filename = 'target_picked.hdf5'

# filepath = path + filename # comment / uncomment this line for batch analysis or single analysis respectively

px = 130


# Parameters that we might want to change

K = 4 # number of NNDs to be calculated
binsize = 0.3 # bin size for NND

nndxlim = 10 # lim of x coord. for NND plot in nm
nndylim = 0.25 # lim of y coord. for NND plot

# create subfolder for results

results_path = os.path.join(os.path.dirname(filepath), 'integrated_results/')
try:
    os.mkdir(results_path)

except OSError:
    pass

# filepath = path + filename 
df_exp = pd.read_hdf(filepath, key = 'locs')

N = df_exp.shape[0] # total number of sugars
print('Total number of sugars:', N)

# N_sim = int(1e6) # this parameter controls the resolution of the simulation
N_sim = int(1e5)

"""
===============================================================================
Simulation of CSR data
===============================================================================
"""

picked_area = io.extract_total_picked_area(filepath.replace('.hdf5', f'.yaml'))
picked_area = np.around(picked_area, 4)

if picked_area is not None:
    print(f"Total Picked Area (μm^2): {np.around(picked_area, 2)}")
else:
    print("The 'Total Picked Area (μm^2)' parameter was not found.")

obs_density = N / picked_area

print('Observed density is:', int(np.around(obs_density, 0)), 'μm^-2')

w = np.sqrt(N_sim/obs_density) 
h = w

# simple CSR model
pos_sim = np.array([np.random.uniform(0, w, N_sim), 
                    np.random.uniform(0, h, N_sim)]).T

pos_sim = pos_sim * 1e3 # in nm

# =============================================================================
# save positions of the simulated molecules as a Picasso compatible file
# =============================================================================

filename_csr = 'simulated_CSR_'

width = w/px
height = h/px

info_csr = {}
info_csr["Generated by"] = "CSR simulation to compare with glycans data"
info_csr["Width"] = width # pixel 
info_csr["Height"] = height # pixel
info_csr["Pixelsize"] = px # in nm

io.save_pos(results_path, filename_csr, width, height, pos_sim/px, [info_csr])

df_sim = pd.read_hdf(results_path + filename_csr + '.hdf5', key = 'locs')

"""
===============================================================================
Analysis pipeline
===============================================================================
"""

clustered_fraction = np.zeros(2) # clustered fraction for 0: exp and 1: CSR
xmax = np.zeros(2) # position of the 1nn peak for  0: exp and 1: CSR

# apply the same analysis to the experimental and the simulated dataset
for i, df in enumerate([df_exp, df_sim]):
    
    print(i)
    # print(df)
    

    """
    ===============================================================================
    Parameters for DBSCAN
    ===============================================================================
    """
    
    epsilon_nm = 10
    epsilon_px = epsilon_nm/px
    minpts = 1
    
    """
    ===============================================================================
    DBSCAN
    ===============================================================================
    """
    # Info to be added to dbscan yaml filename
    # info_db = {
    #     'Generated by': 'DBSCAN',
    #     'epsilon': epsilon_px,
    #     'minpts': minpts
    #     }
        
    db_clusters = dbscan.dbscan_f(df, epsilon_px, minpts)
    
    db_clusters_df = pd.DataFrame.from_records(db_clusters) # convert to an aux dataframe to add the "nsugars in cluster" column
    
    db_clusters_df['nsugars_in_cluster'] =  db_clusters_df.shape[0] * [-1] # -1 means identity not assigned

    """
    ===============================================================================
    Analysis on each cluster (sugar counts, max distance, etc)
    ===============================================================================
    """
    
    nclusters = db_clusters['group'].max()
    
    cluster_size = []
    maxdist_list = []
    
    for j in range(nclusters):
        
        cluster = db_clusters[db_clusters['group'] == j]
        
        sugars = np.array([cluster['x'],
                           cluster['y']])
                    
        nsugars = (sugars.shape[1])
        
        cluster_size.append(nsugars)
        
        db_clusters_df.loc[db_clusters_df['group'] == j, 'nsugars_in_cluster'] = nsugars # add the "nsugars in cluster" column
        
        if nsugars > 2:
        
            pairwise_distances = pdist(sugars.T)
            maxdist = np.max(pairwise_distances)
            
            maxdist_list.append(maxdist)
            
    if i == 0:
    
        info = io.load_info(filepath.replace('.hdf5', f'.yaml'))
        
        db_clusters = db_clusters_df.to_records(index=False) # convert back to recarray
        
        # save locs in dbscan cluster with colorcoding = protein ID
        dbscan_filename = '%s_dbscan_%s_%d.hdf5' % (filename.replace('.hdf5', ''), str(epsilon_nm), minpts)
        io.save_locs(results_path + dbscan_filename, db_clusters, info)
    
    """
    ===============================================================================
    0. Sugar count per cluster
    ===============================================================================
    """
    
    nbins_0 = np.arange(0, 15, 1)
    
    counts0, binedges0 = np.histogram(cluster_size, bins=nbins_0, density=True)
    
    # bin_centers0 = (bin_edges0[:-1] + bin_edges0[1:])/2
    
    if i == 0:
        fig_0, ax_0 = plt.subplots(figsize=(6,5))
    
        ax_0.bar(binedges0[:-1], counts0, edgecolor='black', linewidth=0.2, 
                    width=0.4, alpha=0.8, color='#2880C4', label='Data')
        
        df_save = pd.DataFrame({'counts': np.append(counts0, np.nan), 
                                'binedges': binedges0}) 
        df_save.to_csv(results_path + 'nsugars.csv', index=False)

        
    elif i == 1:
        
        ax_0.scatter(binedges0[2:-1], counts0[2:], edgecolor='black', 
                     linewidth=1, facecolor='None', label='CSR')
        
        df_save = pd.DataFrame({'counts': np.append(counts0, np.nan), 
                                'binedges': binedges0}) 
        df_save.to_csv(results_path + 'nsugars_CSR.csv', index=False)
        
        ax_0.set_xlabel('Sugars per cluster')
        ax_0.set_ylabel('Counts')
        
        ax_0.set_xlim(1.5, 15)
        ax_0.set_ylim(0, counts0[2:].max() * 1.5)
        
        ax_0.legend()
        
    clustered_fraction[i] = counts0[2:].sum()
    
    print('Fraction of clustered sugars:', np.around(counts0[2:].sum(),2))
    
    
    """
    ===============================================================================
    1. Analysis of max dist within a cluster
    ===============================================================================
    """
    
    maxdist_list = np.array(maxdist_list) * px
    
    nbins_1 = np.arange(0, 100, 1)
    
    if i == 0:
        
        fig_1, ax_1 = plt.subplots(figsize=(5,5))
    
        counts1, binedges1, _ = ax_1.hist(maxdist_list, bins=nbins_1, density=True, alpha=0.5, 
                                          label='Data')
        
        df_save = pd.DataFrame({'counts': np.append(counts1, np.nan), 
                                'binedges': binedges1}) 
        df_save.to_csv(results_path + 'maxdist.csv', index=False)
        
    elif i == 1:
        
        counts1, bin_edges1 = np.histogram(maxdist_list, bins=nbins_1, 
                                           density=True)
        
        df_save = pd.DataFrame({'counts': np.append(counts1, np.nan), 
                                'binedges': binedges1}) 
        df_save.to_csv(results_path + 'maxdist_CSR.csv', index=False)
        
        ax_1.plot(bin_edges1[:-1], counts1[:], color='black', 
                  linewidth=1, label='CSR')
        
        ax_1.set_xlabel('Max. distance within a cluster')
        ax_1.set_ylabel('Counts')
        
        ax_1.legend()

    
    
    """
    ===============================================================================
    2. NND calculation for experimental data
    ===============================================================================
    """
    
    # parameters for the NND analysis
    # binsize = 1
    maxdist = 1000
    
    fsize = (10, 10)
    msize = 50
    
    x_0 = df.x * px # in nm
    y_0 = df.y * px # in nm
    
    pos_exp_0 = np.array([x_0, y_0]).T
    
    # find nearest neighbours from biomolecule 0 to biomolecule 1 (exp data), in this case is the same biomolecule (e.g. ManNAz, GalNAz)
    nbrs = NearestNeighbors(n_neighbors=5).fit(pos_exp_0) 
    _distances_exp, _indices_exp = nbrs.kneighbors(pos_exp_0) # get distances and indices
    
    distances_exp = {}
    freqs_exp = {}
    bin_centers_exp = {}
    
    bins = np.arange(0, maxdist, binsize)

    for j in range(K):
        
        key = str(j+1) + 'nn'
        distances_exp[key] = _distances_exp[:, j+1] # get the first neighbour distances (0 for hetero, 1 for homo in the second coord)
    
        freqs_exp[key], binedges = np.histogram(distances_exp[key], bins=bins, 
                                                density=True)
            
        bin_centers_exp[key] = (binedges[:-1] + binedges[1:])/2
        
        if j == 0: # apply KDE to the 1nn
            
            data = distances_exp[key]
            xxkde = np.linspace(0, 200, 2000)
                            
            kde = gaussian_kde(data, bw_method=0.1)  # bw_method sets the bandwidth
            
    if i == 0:
        
        for j in range(K):
            
            key = str(j+1) + 'nn'
            freqs_exp[key] = np.append(freqs_exp[key], np.nan)
        
        df_save = pd.DataFrame(freqs_exp)
        df_save['binedges'] = binedges
        
        df_save.to_csv(results_path + 'nnds.csv', index=False)
        
        df_save_distances = pd.DataFrame(distances_exp)
        df_save_distances.to_csv(results_path + 'distances.csv', index=False)
        
    if i == 1:
        
        for j in range(K):
            
            key = str(j+1) + 'nn'
            freqs_exp[key] = np.append(freqs_exp[key], np.nan)
        
        df_save = pd.DataFrame(freqs_exp)
        df_save['binedges'] = binedges
        
        df_save.to_csv(results_path + 'nnds_csr.csv', index=False)
        
        df_save_distances = pd.DataFrame(distances_exp)
        df_save_distances.to_csv(results_path + 'distances_csr.csv', index=False)
        

    # =============================================================================
    # Plot of experimental NNDs
    # =============================================================================
    
    colors = ['#2880C4', '#97D8C4', '#F4B942', '#363636']
    
    if i == 0:
    
        fig_2, ax_2 = plt.subplots(figsize=(5,5))
        
        # ax_2.set_xlim(0, 100)
        # ax_2.set_ylim(0, np.nanmax(freqs_exp['1nn']) * 1.5)
        
        ax_2.set_xlim(0, nndxlim)
        ax_2.set_ylim(0, nndylim)
        
        ax_2.set_xlabel('K-th NND (nm)')
        ax_2.set_ylabel('Counts')
        
        ax_2.set_box_aspect(1)
        
        
        for j in range(K):
            
            key = str(j+1) + 'nn'
        
            counts, bin_edges, _ = ax_2.hist(distances_exp[key], bins=bins, 
                                             edgecolor='black', linewidth=0.1, 
                                             alpha=0.5, density=True, 
                                             color=colors[j], label=key)
            
            if j == 0:
                
                ax_2.plot(xxkde, kde(xxkde), 'k--')
                xmax[i] = xxkde[np.argmax(kde(xxkde))]
        
                        
    elif i == 1:
        
        for j in range(K):
            
            key = str(j+1) + 'nn'
        
            counts, binedges = np.histogram(distances_exp[key], 
                                            bins=bins, density=True)
                                                
            ax_2.plot(binedges[:-1], counts, color=colors[j])
            
            if j == 0:
                
                ax_2.plot(xxkde, kde(xxkde), 'k--')
                xmax[i] = xxkde[np.argmax(kde(xxkde))]
            
nnd_peak_shift = (1 - (xmax[1] - xmax[0]) / xmax[1]) * 100

print('NND peak shift (rel. to CSR) is', np.around(nnd_peak_shift, 1), '%')

clustered_fraction_increase = (clustered_fraction[0]/clustered_fraction[1] - 1) * 100

print('Clustered fraction increase (rel. to CSR) is', np.around(clustered_fraction_increase, 1), '%')

            
fig_0.savefig(results_path + 'nsugars.pdf', format='pdf')
fig_1.savefig(results_path + 'maxdists.pdf', format='pdf')
fig_2.savefig(results_path + 'nnds.pdf', format='pdf')


final_numbers = {
    "obs_density (μm^-2)": [np.around(obs_density, 2)],
    "clustered_fraction (rel. increase in %)": [np.around(clustered_fraction_increase, 2)],
    "nnd_peak_shift (rel. increase in %) ": [np.around(nnd_peak_shift, 2)]
}

# Create a DataFrame and save to CSV
fn = pd.DataFrame(final_numbers)
fn.to_csv(results_path + 'final_numbers.csv', index=False)
            
            