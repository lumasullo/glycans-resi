#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 15:21:59 2024

@author: masullo
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set main directory
main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/ManNAz/ManNAz_combined'
basename = os.path.basename(main_dir)
dataset = basename.split('_combined')[0]
print(f"Dataset: {dataset}")

# Define the distance metrics to analyze
distance_metrics = ['1nn', '2nn', '3nn', '4nn']

# Generate a random subset of 40 numbered folders (1 to 40)

# import random
# subset = random.sample(range(1, 41), 40)

subset = range(1, 40)

string_list = [f"{num}" for num in subset]
print("Selected folders:", string_list)

# Initialize dictionaries to store data arrays for each metric
all_distances = {metric: [] for metric in distance_metrics}
all_distances_csr = {metric: [] for metric in distance_metrics}

# Loop through each numbered folder in the main directory
for folder_name in sorted(string_list):
    folder_path = os.path.join(main_dir, folder_name)
    print(f"Processing folder: {folder_name}")

    if os.path.isdir(folder_path):
        csv_path_distances = os.path.join(folder_path, 'integrated_results/distances.csv')
        csv_path_distances_csr = os.path.join(folder_path, 'integrated_results/distances_csr.csv')

        # Process distances.csv
        if os.path.exists(csv_path_distances):
            df_distances = pd.read_csv(csv_path_distances)
            for metric in distance_metrics:
                if metric in df_distances.columns:
                    all_distances[metric].append(df_distances[metric].values)
                else:
                    print(f"'{metric}' not found in {csv_path_distances}")

        # Process distances_csr.csv
        if os.path.exists(csv_path_distances_csr):
            df_distances_csr = pd.read_csv(csv_path_distances_csr)
            for metric in distance_metrics:
                if metric in df_distances_csr.columns:
                    all_distances_csr[metric].append(df_distances_csr[metric].values)
                else:
                    print(f"'{metric}' not found in {csv_path_distances_csr}")

# Concatenate data arrays for each metric
combined_distances = {metric: np.concatenate(all_distances[metric]) if all_distances[metric] else np.array([]) 
                      for metric in distance_metrics}
combined_distances_csr = {metric: np.concatenate(all_distances_csr[metric]) if all_distances_csr[metric] else np.array([]) 
                          for metric in distance_metrics}

# Plot settings
colors = ['#2880C4', '#F4B942', '#D9534F', '#5CB85C']  # Colors for NN1, NN2, NN3, NN4
labels = ['1NN', '2NN', '3NN', '4NN']
binsize = 0.3
maxdist = 100
nndxlim = 30
nndylim = 0.09
bins = np.arange(0, maxdist, binsize)

# Merged Histograms for Data and CSR (Line Plots for CSR)
fig, ax = plt.subplots(figsize=(10, 6))

ax.set_xlim(0, nndxlim)
ax.set_ylim(0, nndylim)
ax.set_xlabel('NND (nm)')
ax.set_ylabel('Frequency')
ax.set_title(f'Merged Histogram - Data and CSR')

# Iterate through each nearest neighbor metric
for idx, metric in enumerate(distance_metrics):
    color = colors[idx]

    # Data Histograms (Bar Plots)
    data_distances = combined_distances[metric]
    if len(data_distances) > 0:
        ax.hist(data_distances, bins=bins, edgecolor='black', linewidth=0.1, alpha=0.5,
                density=True, color=color, label=f'{labels[idx]} data')

    # CSR Histograms (Line Plots)
    data_distances_csr = combined_distances_csr[metric]
    if len(data_distances_csr) > 0:
        # Calculate the bin centers for CSR line plot
        bin_centers = (bins[:-1] + bins[1:]) / 2
        counts, _ = np.histogram(data_distances_csr, bins=bins, density=True)
        ax.plot(bin_centers, counts, label=f'{labels[idx]} CSR', color=color, linewidth=2)

# Add legends
ax.legend()

# Show plots
plt.tight_layout()
plt.show()

fig.savefig(main_dir + '/output.pdf', format='pdf', bbox_inches='tight', pad_inches=0.1, transparent=True)


