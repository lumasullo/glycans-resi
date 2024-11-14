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

# Histogram and KDE for combined_1nn_array_distances
freqs_exp_distances, binedges_distances = np.histogram(data_distances, bins=bins, density=True)
bin_centers_exp_distances = (binedges_distances[:-1] + binedges_distances[1:]) / 2

# Histogram and KDE for combined_1nn_array_distances_csr
freqs_exp_distances_csr, binedges_distances_csr = np.histogram(data_distances_csr, bins=bins, density=True)
bin_centers_exp_distances_csr = (binedges_distances_csr[:-1] + binedges_distances_csr[1:]) / 2

# KDE for both datasets
xxkde = np.linspace(0, 200, 2000)
kde_distances = gaussian_kde(data_distances, bw_method=0.1)
kde_distances_csr = gaussian_kde(data_distances_csr, bw_method=0.1)

# Plotting
fig_0, ax_0 = plt.subplots(figsize=(5, 5))

# General plot settings
ax_0.set_xlim(0, nndxlim)
ax_0.set_ylim(0, nndylim)
ax_0.set_xlabel('1st NND (nm)')
ax_0.set_ylabel('Freq.')
ax_0.set_box_aspect(1)

# Histogram and KDE plot for data_distances
counts, bin_edges, _ = ax_0.hist(data_distances, bins=bins, edgecolor='black', linewidth=0.1,
                                 alpha=0.5, density=True, color=colors[0], label='data')
ax_0.plot(xxkde, kde_distances(xxkde), linestyle='--', color=colors[0])

# Histogram and KDE plot for data_distances_csr
counts_csr, bin_edges_csr, _ = ax_0.hist(data_distances_csr, bins=bins, edgecolor='black', linewidth=0.1,
                                         alpha=0.2, density=True, color=colors[1], label='CSR')
ax_0.plot(xxkde, kde_distances_csr(xxkde), linestyle='--', color=colors[1])

# Adding legend for clarity
ax_0.legend()

# Display the plot
plt.show()


"""
===============================================================================
3. NND analysis with mean and std error of the mean for the histogram
===============================================================================
"""

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


# Define the range of interest for integration (0 to 10 nm)
range_min, range_max = 0, 10

# Find the indices of bin_centers that fall within this range
indices_within_range = np.where((bin_centers >= range_min) & (bin_centers <= range_max))[0]

# Calculate the integrated area between the mean curves within the 0-10 nm range
area_between_means = np.trapz(np.abs(mean_heights_distances[indices_within_range] - mean_heights_distances_csr[indices_within_range]), 
                              bin_centers[indices_within_range])

# Calculate the uncertainty by integrating the area between the std deviation curves within the 0-10 nm range
area_uncertainty = np.trapz(std_heights_distances[indices_within_range] + std_heights_distances_csr[indices_within_range], 
                            bin_centers[indices_within_range])

# Plotting the mean curves with shaded area between them in the range of 0 to 10 nm
fig_2, ax_2 = plt.subplots()

ax_2.set_title('1st NND histogram and area - ' + dataset)

# Plot for all_1nn_arrays_distances
ax_2.plot(bin_centers, mean_heights_distances, linestyle='-', color=colors[0], label='Mean Height (distances)')
ax_2.plot(bin_centers, mean_heights_distances + std_heights_distances, linestyle='--', color=colors[0], label='+1 Std Dev (distances)')
ax_2.plot(bin_centers, mean_heights_distances - std_heights_distances, linestyle='--', color=colors[0], label='-1 Std Dev (distances)')

# Plot for all_1nn_arrays_distances_csr
ax_2.plot(bin_centers, mean_heights_distances_csr, linestyle='-', color=colors[1], label='Mean Height (distances_csr)')
ax_2.plot(bin_centers, mean_heights_distances_csr + std_heights_distances_csr, linestyle='--', color=colors[1], label='+1 Std Dev (distances_csr)')
ax_2.plot(bin_centers, mean_heights_distances_csr - std_heights_distances_csr, linestyle='--', color=colors[1], label='-1 Std Dev (distances_csr)')

# Add light gray shading between the two mean curves within the 0-10 nm range
ax_2.fill_between(bin_centers[indices_within_range],
                  mean_heights_distances[indices_within_range],
                  mean_heights_distances_csr[indices_within_range],
                  color='lightgray', alpha=0.4)

# General plot settings
ax_2.set_xlim(0, nndxlim)
ax_2.set_ylim(0, nndylim)
ax_2.set_xlabel('1st NND (nm)')
ax_2.set_ylabel('Freq.')
ax_2.set_box_aspect(1)

# Plotting the integrated area with uncertainty as a bar plot
fig_3, ax_3 = plt.subplots()

ax_3.bar(['Integrated Area (0-10 nm)'], [area_between_means], yerr=[area_uncertainty], color='gray', alpha=0.7, capsize=10, width=0.3)
ax_3.set_ylabel('Integrated Area')
ax_3.set_title('Integrated Area between Mean Curves (0-10 nm) - ' + dataset)
ax_3.set_xlim(-0.5, 0.5)
ax_3.set_ylim(0, 0.27)
ax_3.set_box_aspect(1)


plt.show()

save_path = main_dir

# Save each figure as a PDF
fig_1.savefig(os.path.join(save_path, 'fig_1.pdf'), format='pdf')
fig_2.savefig(os.path.join(save_path, 'fig_2.pdf'), format='pdf')
fig_3.savefig(os.path.join(save_path, 'fig_3.pdf'), format='pdf')



