#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:57:35 2024

@author: masullo
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.close()

# Directories for ManNAz and GalNAz data
main_dirs = [
    '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/ManNAz/ManNAz_combined',
    '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/GalNAz/GalNAz_combined'
]

# Initialize lists to store counts data for each file (for both datasets)
all_counts_nsugars_manna = []
all_counts_nsugars_csr_manna = []
all_counts_nsugars_galna = []
all_counts_nsugars_csr_galna = []

# Loop through each directory
for main_dir in main_dirs:
    # Initialize temporary lists for the current directory
    all_counts_nsugars = []
    all_counts_nsugars_csr = []

    # Loop through each numbered folder in the current main directory
    for folder_name in sorted(os.listdir(main_dir)):
        folder_path = os.path.join(main_dir, folder_name)
        
        # Check if the folder path is a directory
        if os.path.isdir(folder_path):
            # Paths for integrated_results/nsugars.csv and integrated_results/nsugars_csr.csv
            csv_path_nsugars = os.path.join(folder_path, 'integrated_results/nsugars.csv')
            csv_path_nsugars_csr = os.path.join(folder_path, 'integrated_results/nsugars_csr.csv')
            
            # Process nsugars.csv
            if os.path.exists(csv_path_nsugars):
                df_nsugars = pd.read_csv(csv_path_nsugars)
                counts_nsugars = df_nsugars.iloc[:, 0].dropna().values  # Remove NaN values
                all_counts_nsugars.append(counts_nsugars)
                # Use bin edges from the first file since they're identical across files
                if 'bin_edges_nsugars' not in locals():
                    bin_edges_nsugars = df_nsugars.iloc[:, 1].dropna().values
            
            # Process nsugars_csr.csv
            if os.path.exists(csv_path_nsugars_csr):
                df_nsugars_csr = pd.read_csv(csv_path_nsugars_csr)
                counts_nsugars_csr = df_nsugars_csr.iloc[:, 0].dropna().values  # Remove NaN values
                all_counts_nsugars_csr.append(counts_nsugars_csr)
                if 'bin_edges_nsugars_csr' not in locals():
                    bin_edges_nsugars_csr = df_nsugars_csr.iloc[:, 1].dropna().values

    # Ensure all arrays have the same length by trimming to the minimum length
    min_len = min(len(arr) for arr in all_counts_nsugars)
    all_counts_nsugars = np.array([arr[:min_len] for arr in all_counts_nsugars])

    min_len_csr = min(len(arr) for arr in all_counts_nsugars_csr)
    all_counts_nsugars_csr = np.array([arr[:min_len_csr] for arr in all_counts_nsugars_csr])

    # Calculate the mean and standard deviation for each bin
    mean_counts_nsugars = np.mean(all_counts_nsugars, axis=0)
    std_counts_nsugars = np.std(all_counts_nsugars, axis=0)

    mean_counts_nsugars_csr = np.mean(all_counts_nsugars_csr, axis=0)
    std_counts_nsugars_csr = np.std(all_counts_nsugars_csr, axis=0)

    # Append results for each dataset (ManNAz and GalNAz)
    if main_dir == '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/ManNAz/ManNAz_combined':
        all_counts_nsugars_manna.append((mean_counts_nsugars, std_counts_nsugars))
        all_counts_nsugars_csr_manna.append((mean_counts_nsugars_csr, std_counts_nsugars_csr))
    else:
        all_counts_nsugars_galna.append((mean_counts_nsugars, std_counts_nsugars))
        all_counts_nsugars_csr_galna.append((mean_counts_nsugars_csr, std_counts_nsugars_csr))

# Plotting
fig, ax = plt.subplots(figsize=(6, 5))

delta_plot = 0.05

# Define colors for easy modification
# colors = {
#     'manna_nsugars': '#404040',  # Dark gray for ManNAz nsugars
#     'manna_nsugars_csr': '#404040',  
#     'galna_nsugars': '#B0B0B0',  
#     'galna_nsugars_csr': '#B0B0B0'  # Light gray for GalNAz CSR
# }

colors = {
    'manna_nsugars': '#FF3C38',  # Dark gray for ManNAz nsugars
    'manna_nsugars_csr': '#FF3C38',  
    'galna_nsugars': '#6C8EAD',  
    'galna_nsugars_csr': '#6C8EAD'  # Light gray for GalNAz CSR
}

K = 1 # change to 41 or the actual number of areas used if s.e.m. is desired

# Plot mean counts with error bars for ManNAz (nsugars)
for i, (mean_counts, std_counts) in enumerate(all_counts_nsugars_manna):
    ax.errorbar(
        bin_edges_nsugars[:-1] + delta_plot, mean_counts, yerr=std_counts/np.sqrt(K), fmt='-o', 
        color=colors['manna_nsugars'], label='ManNAz data', capsize=3
    )

# Plot solid markers with transparency for ManNAz (nsugars_csr)
for i, (mean_counts, std_counts) in enumerate(all_counts_nsugars_csr_manna):
    ax.errorbar(
        bin_edges_nsugars_csr[:-1] + delta_plot, mean_counts, yerr=std_counts/np.sqrt(K), fmt='--^', 
        color=colors['manna_nsugars_csr'], alpha=0.7, label='ManNAz CSR', capsize=3
    )

# Plot mean counts with error bars for GalNAz (nsugars)
for i, (mean_counts, std_counts) in enumerate(all_counts_nsugars_galna):
    ax.errorbar(
        bin_edges_nsugars[:-1] - delta_plot, mean_counts, yerr=std_counts/np.sqrt(K), fmt='-o', 
        color=colors['galna_nsugars'], label='GalNAz data', capsize=3
    )

# Plot solid markers with transparency for GalNAz (nsugars_csr)
for i, (mean_counts, std_counts) in enumerate(all_counts_nsugars_csr_galna):
    ax.errorbar(
        bin_edges_nsugars_csr[:-1] - delta_plot, mean_counts, yerr=std_counts/np.sqrt(K), fmt='--^', 
        color=colors['galna_nsugars_csr'], alpha=0.7, label='GalNAz CSR', capsize=3
    )

# Set y-axis to logarithmic scale
ax.set_yscale('log')

# Plot settings
ax.set_xlim(1.8, 5.2)
ax.set_ylim(1e-4, 1)
ax.set_xlabel('N sugars per cluster')
ax.set_ylabel('Counts')

# Set x-axis ticks at 2, 3, 4, 5 and format them as integers
ax.set_xticks([2, 3, 4, 5])
ax.set_xticklabels([2, 3, 4, 5], fontsize=12)
# ax.set_title('Mean Counts with Standard Deviation')

# Combine the legends for the labels
ax.legend()

plt.show()

integrated_csr_manna = all_counts_nsugars_csr_manna[0][0][2:-2].sum()
integrated_manna = all_counts_nsugars_manna[0][0][2:-2].sum()

delta_manna = integrated_manna/integrated_csr_manna

integrated_csr_galna = all_counts_nsugars_csr_galna[0][0][2:-2].sum()
integrated_galna = all_counts_nsugars_galna[0][0][2:-2].sum()

delta_galna = integrated_galna/integrated_csr_galna

