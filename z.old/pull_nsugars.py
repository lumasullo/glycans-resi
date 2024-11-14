#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 13:44:52 2024

@author: masullo
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.close()

# Set the main directory containing the numbered folders
# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/GalNAz/GalNAz_combined'
main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/ManNAz/ManNAz_combined'


# Initialize lists to store counts data for each file
all_counts_nsugars = []
all_counts_nsugars_csr = []

# Loop through each numbered folder in the main directory
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

# Calculate bin centers, ensuring they match the length of the counts arrays
# bin_centers_nsugars = (bin_edges_nsugars[:len(mean_counts_nsugars)] + bin_edges_nsugars[1:len(mean_counts_nsugars) + 1]) / 2
# bin_centers_nsugars_csr = (bin_edges_nsugars_csr[:len(mean_counts_nsugars_csr)] + bin_edges_nsugars_csr[1:len(mean_counts_nsugars_csr) + 1]) / 2

# Plotting
fig, ax = plt.subplots(figsize=(6, 5))

# Plot mean counts with error bars for nsugars (dots only, with caps)
ax.errorbar(
    bin_edges_nsugars[:-1], mean_counts_nsugars, yerr=std_counts_nsugars, fmt='o',
    color='#2880C4', label='nsugars.csv', capsize=3
)

# Plot mean counts with error bars for nsugars_csr (dots only, with caps)
ax.errorbar(
   bin_edges_nsugars_csr[:-1], mean_counts_nsugars_csr, yerr=std_counts_nsugars_csr, fmt='o',
    color='#F4B942', label='nsugars_csr.csv', capsize=3
)

# Set y-axis to logarithmic scale
ax.set_yscale('log')

ax.set_xlim(1.8, 5.2)
ax.set_ylim(1e-4, 1)

# Plot settings
ax.set_xlabel('Distance (nm)')
ax.set_ylabel('Counts')
ax.set_title('Mean Counts with Standard Deviation')
ax.legend()

plt.show()
