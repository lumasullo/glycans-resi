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

main_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/GalNAz/GalNAz_combined'
# main_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz/ManNAz_combined'


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

# yellow, #F0E442, 240, 228, 66 ; blue, #0072B2, 0, 114, 178 ; vermilion, #D55E00, 213, 94, 0 ; reddish purple, #CC79A7, 204, 121, 167


# colors = ['#2880C4', '#F4B942', '#D9534F', '#5CB85C']  # Colors for NN1, NN2, NN3, NN4
# colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442']
colors = ['#2880C4', '#F4B942', '#D9534F', '#009E73']  # Colors for NN1, NN2, NN3, NN4

labels = ['1NN', '2NN', '3NN', '4NN']
binsize = 0.25
maxdist = 100
nndxlim = 30
nndylim = 0.08
bins = np.arange(0, maxdist, binsize)

# Merged Histograms for Data and CSR (Line Plots for CSR)
fig, ax = plt.subplots(figsize=(9, 4))

ax.set_xlim(0, nndxlim)
ax.set_ylim(0, nndylim)
ax.set_xlabel('NND (nm)')
ax.set_ylabel('Frequency')
ax.set_title(f'Merged Histogram - Data and CSR')
ax.tick_params(axis='both', which='both', direction='in')

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
# ax.legend()

# Show plots
plt.tight_layout()
plt.show()

fig.savefig(main_dir + '/kth_nnds.pdf', format='pdf', bbox_inches='tight', pad_inches=0.1, transparent=True)

# === Save combined data (raw) ===

data_save_path = os.path.join(main_dir, 'combined_distances_data_only.csv')
rows_data = []

for metric in distance_metrics:
    for val in combined_distances[metric]:
        rows_data.append({'metric': metric, 'value': val})

df_data = pd.DataFrame(rows_data)
df_data.to_csv(data_save_path, index=False)
print(f"Saved raw data to {data_save_path}")

# === Save CSR line plot data (bin centers + density) ===

csr_save_path = os.path.join(main_dir, 'csr_line_data.csv')
rows_csr = []

for metric in distance_metrics:
    csr_vals = combined_distances_csr[metric]
    if len(csr_vals) > 0:
        counts, _ = np.histogram(csr_vals, bins=bins, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        for x, y in zip(bin_centers, counts):
            rows_csr.append({'metric': metric, 'bin_center': x, 'density': y})

df_csr = pd.DataFrame(rows_csr)
df_csr.to_csv(csr_save_path, index=False)
print(f"Saved CSR line plot data to {csr_save_path}")

# === Reload and replot from saved CSVs ===

df_data_loaded = pd.read_csv(data_save_path)
df_csr_loaded = pd.read_csv(csr_save_path)

fig, ax = plt.subplots(figsize=(9, 4))
ax.set_xlim(0, nndxlim)
ax.set_ylim(0, nndylim)
ax.set_xlabel('NND (nm)')
ax.set_ylabel('Frequency')
ax.set_title('retrieved from saved data and CSR line')
ax.tick_params(axis='both', which='both', direction='in')

for idx, metric in enumerate(distance_metrics):
    color = colors[idx]

    # Plot data histogram
    vals = df_data_loaded[df_data_loaded['metric'] == metric]['value'].values
    if len(vals) > 0:
        ax.hist(vals, bins=bins, edgecolor='black', linewidth=0.1, alpha=0.5,
                density=True, color=color, label=f'{labels[idx]} data')

    # Plot CSR line
    df_csr_metric = df_csr_loaded[df_csr_loaded['metric'] == metric]
    if not df_csr_metric.empty:
        ax.plot(df_csr_metric['bin_center'], df_csr_metric['density'],
                label=f'{labels[idx]} CSR', color=color, linewidth=2)

plt.tight_layout()
plt.show()



