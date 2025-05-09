#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 10:31:31 2024

@author: masullo
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Set the main directory containing the numbered folders
# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/Glycosylated spherical domains HMECs/Spherical clusters/240618HMECmannaz/'
# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/Homogenous areas/GalNAz/240618_GalNAz'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/Homogenous areas/ManNAz/240618_ManNAz'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/Glycosylated spherical domains HMECs/Spherical clusters/240617galnaz'
# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/Glycosylated spherical domains HMECs/Spherical clusters/240618galnaz'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/High density nanoclusters/ManNAz/240618_ManNAz'
# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/High density nanoclusters/GalNAz/240618_GalNAz'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/MCF10A_and_MCF10AT/Spherical MCF10As/240819_MCF10A'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/MCF10A_and_MCF10AT/Spherical MCF10As/240820_MCF10AT'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/spherical nanodomains/20240618_ManNAz'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/MCF10A_and_MCF10AT/Homogenous  ROIs/240820_MCF10AT'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/MCF10A_and_MCF10AT/High density domains/240820_MCF10AT'

main_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz/ManNAz_combined'

# main_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/GalNAz/GalNAz_combined'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/MCF10As Homogenous  ROIs/MCF10A_combined'

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/MCF10As Homogenous  ROIs/MCF10AT_combined'

# Extract the part before "_combined"
basename = os.path.basename(main_dir)  # Get the last part of the path
dataset = basename.split('_combined')[0]

print(dataset)


import random

# Generate a random subset of 10 elements from numbers 1 to 40
subset = random.sample(range(1, 41), 40)

# Create a list of strings, one for each element in the subset
string_list = [f"{num}" for num in subset]

# Print the results
print("Random subset:", subset)
print("List of strings:", string_list)


# Initialize empty lists to store all "1nn" arrays for each file
all_1nn_arrays_distances = []
all_1nn_arrays_distances_csr = []

# Loop through each numbered folder in the main directory
# for folder_name in sorted(os.listdir(main_dir)):
    
for folder_name in sorted(string_list ):
    
    folder_path = os.path.join(main_dir, folder_name)
    
    print(folder_name)
    
    # Check if the folder path is a directory
    if os.path.isdir(folder_path):
        # Paths for distances.csv and distances_csr.csv
        csv_path_distances = os.path.join(folder_path, 'integrated_results/distances.csv')
        csv_path_distances_csr = os.path.join(folder_path, 'integrated_results/distances_csr.csv')
        
        # Process distances.csv
        if os.path.exists(csv_path_distances):
            df_distances = pd.read_csv(csv_path_distances)
            if '1nn' in df_distances.columns:
                all_1nn_arrays_distances.append(df_distances['1nn'].values)
            else:
                print(f"'1nn' not found in {csv_path_distances}")
        else:
            print(f"{csv_path_distances} does not exist.")
        
        # Process distances_csr.csv
        if os.path.exists(csv_path_distances_csr):
            df_distances_csr = pd.read_csv(csv_path_distances_csr)
            if '1nn' in df_distances_csr.columns:
                all_1nn_arrays_distances_csr.append(df_distances_csr['1nn'].values)
            else:
                print(f"'1nn' not found in {csv_path_distances_csr}")
        else:
            print(f"{csv_path_distances_csr} does not exist.")

# Concatenate all "1nn" arrays into single arrays for each file
combined_1nn_array_distances = np.concatenate(all_1nn_arrays_distances) if all_1nn_arrays_distances else np.array([])
combined_1nn_array_distances_csr = np.concatenate(all_1nn_arrays_distances_csr) if all_1nn_arrays_distances_csr else np.array([])

# Print the results for verification
print("Combined 1nn array from distances.csv:", combined_1nn_array_distances)
print("Combined 1nn array from distances_csr.csv:", combined_1nn_array_distances_csr)


"""
===============================================================================
2. NND analysis of the pulled 1st NND data
===============================================================================
"""

colors = ['#2880C4', '#404040']  # Adjusted to differentiate the two arrays
binsize = 0.3
maxdist = 100
nndxlim = 15 
nndylim = 0.09
bins = np.arange(0, maxdist, binsize)

# Define data arrays
data_distances = combined_1nn_array_distances
data_distances_csr = combined_1nn_array_distances_csr

# Function to calculate mean and std heights for a list of arrays
def calculate_mean_std_histograms(data_list, bins, N):
    histograms = []
    for arr in data_list:
        hist, _ = np.histogram(arr, bins=bins, density=True)
        histograms.append(hist)
    histograms = np.array(histograms)
    mean_heights = np.mean(histograms, axis=0)
    std_heights = np.std(histograms, axis=0) / np.sqrt(N) # std of the mean
    return mean_heights, std_heights

N = len(all_1nn_arrays_distances) # get number of areas analyzed

# Calculate mean and std for each list
mean_heights_distances, std_heights_distances = calculate_mean_std_histograms(all_1nn_arrays_distances, bins, N)
mean_heights_distances_csr, std_heights_distances_csr = calculate_mean_std_histograms(all_1nn_arrays_distances_csr, bins, N)

# Plotting
bin_centers = (bins[:-1] + bins[1:]) / 2  # Calculate bin centers for plotting

fig_1, ax_1 = plt.subplots()

# Plot for all_1nn_arrays_distances
ax_1.plot(bin_centers, mean_heights_distances, linestyle='-', color=colors[0], label='Mean Height (distances)')
ax_1.plot(bin_centers, mean_heights_distances + std_heights_distances, linestyle='--', color=colors[0], label='+1 Std Dev (distances)')
ax_1.plot(bin_centers, mean_heights_distances - std_heights_distances, linestyle='--', color=colors[0], label='-1 Std Dev (distances)')
ax_1.fill_between(bin_centers, mean_heights_distances + std_heights_distances, mean_heights_distances - std_heights_distances, 
                  color=colors[0], alpha=0.3)

# Plot for all_1nn_arrays_distances_csr
ax_1.plot(bin_centers, mean_heights_distances_csr, linestyle='-', color=colors[1], label='Mean Height (distances_csr)')
ax_1.plot(bin_centers, mean_heights_distances_csr + std_heights_distances_csr, linestyle='--', color=colors[1], label='+1 Std Dev (distances_csr)')
ax_1.plot(bin_centers, mean_heights_distances_csr - std_heights_distances_csr, linestyle='--', color=colors[1], label='-1 Std Dev (distances_csr)')
ax_1.fill_between(bin_centers, mean_heights_distances_csr + std_heights_distances_csr, mean_heights_distances_csr - std_heights_distances_csr, 
                  color=colors[1], alpha=0.2)

# General plot settings
ax_1.set_xlim(0, nndxlim)
ax_1.set_ylim(0, nndylim)
ax_1.set_xlabel('K-th NND (nm)')
ax_1.set_ylabel('Freq.')
ax_1.set_box_aspect(1)

ax_1.set_title('1st NND histogram - ' + dataset)


plt.show()

save_path = main_dir

# Save each figure as a PDF
fig_1.savefig(os.path.join(save_path, 'fig_1.pdf'), format='pdf')
fig_2.savefig(os.path.join(save_path, 'fig_2.pdf'), format='pdf')
fig_3.savefig(os.path.join(save_path, 'fig_3.pdf'), format='pdf')


