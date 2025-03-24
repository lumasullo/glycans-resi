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
    '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz/ManNAz_combined',
    '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/GalNAz/GalNAz_combined'
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
    if main_dir == '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz/ManNAz_combined':
        all_counts_nsugars_manna.append((mean_counts_nsugars, std_counts_nsugars))
        all_counts_nsugars_csr_manna.append((mean_counts_nsugars_csr, std_counts_nsugars_csr))
    else:
        all_counts_nsugars_galna.append((mean_counts_nsugars, std_counts_nsugars))
        all_counts_nsugars_csr_galna.append((mean_counts_nsugars_csr, std_counts_nsugars_csr))


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

K = 40 # change to 41 or the actual number of areas used if s.e.m. is desired

# Plotting
fig, ax = plt.subplots(figsize=(6, 5))

# Define bar width
bar_width = 0.4

# Positions for each group of bars
n_values = [2, 3, 4, 5]  # Bin edges for N sugars per cluster
x = np.arange(len(n_values))  # Positions for bar groups

# Define alternative colors for CSR bars
colors = {
    'manna_nsugars': '#FF3C38',       # Original red for ManNAz
    'manna_nsugars_csr': '#FF9999',  # Lighter red for ManNAz CSR
    'galna_nsugars': '#6C8EAD',      # Original blue for GalNAz
    'galna_nsugars_csr': '#99CCFF'   # Lighter blue for GalNAz CSR
}

# Plot bars for ManNAz
mean_manna, std_manna = all_counts_nsugars_manna[0]
ax.bar(
    x - bar_width / 2, 
    mean_manna[2:6], 
    bar_width, 
    yerr=std_manna[2:6] / np.sqrt(K), 
    label='ManNAz', 
    color=colors['manna_nsugars'], 
    capsize=3,
    bottom = 1e-4
)

# Plot CSR bars for ManNAz
mean_manna_csr, std_manna_csr = all_counts_nsugars_csr_manna[0]
ax.bar(
    x - bar_width / 2, 
    mean_manna_csr[2:6], 
    bar_width, 
    yerr=std_manna_csr[2:6] / np.sqrt(K), 
    label='ManNAz CSR', 
    color=colors['manna_nsugars_csr'], 
    capsize=3,
    bottom = 1e-4
)

# Plot bars for GalNAz
mean_galna, std_galna = all_counts_nsugars_galna[0]
ax.bar(
    x + bar_width / 2, 
    mean_galna[2:6], 
    bar_width, 
    yerr=std_galna[2:6] / np.sqrt(K), 
    label='GalNAz', 
    color=colors['galna_nsugars'], 
    capsize=3,
    bottom = 1e-4
)

# Plot CSR bars for GalNAz
mean_galna_csr, std_galna_csr = all_counts_nsugars_csr_galna[0]
ax.bar(
    x + bar_width / 2, 
    mean_galna_csr[2:6], 
    bar_width, 
    yerr=std_galna_csr[2:6] / np.sqrt(K), 
    label='GalNAz CSR', 
    color=colors['galna_nsugars_csr'], 
    capsize=3,
    bottom = 1e-4
)


# Set y-axis to logarithmic scale
ax.set_yscale('log')

# Plot settings
ax.set_xlim(-0.5, len(x) - 0.5)  # Adjust x-axis limits
ax.set_ylim(1e-4, 1)  # Adjust y-axis limits based on your data
ax.set_xlabel('N sugars per cluster')
ax.set_ylabel('Counts')

# Set x-axis ticks and labels
ax.set_xticks(x)
ax.set_xticklabels(n_values, fontsize=12)

# Combine the legends for the labels
ax.legend()

plt.tight_layout()
plt.show()



