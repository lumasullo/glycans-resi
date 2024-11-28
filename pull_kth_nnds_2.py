#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:04:53 2024

@author: masullo
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from sklearn.utils import resample

# Set main directory
main_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz/ManNAz_combined'
basename = os.path.basename(main_dir)
dataset = basename.split('_combined')[0]
print(f"Dataset: {dataset}")

# Define the distance metrics to analyze
distance_metrics = ['1nn', '2nn', '3nn', '4nn']

# Generate a random subset of 40 numbered folders (1 to 40)
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

        # Process distances_csr.csv
        if os.path.exists(csv_path_distances_csr):
            df_distances_csr = pd.read_csv(csv_path_distances_csr)
            for metric in distance_metrics:
                if metric in df_distances_csr.columns:
                    all_distances_csr[metric].append(df_distances_csr[metric].values)

# Concatenate data arrays for each metric
combined_distances = {metric: np.concatenate(all_distances[metric]) if all_distances[metric] else np.array([]) 
                      for metric in distance_metrics}
combined_distances_csr = {metric: np.concatenate(all_distances_csr[metric]) if all_distances_csr[metric] else np.array([]) 
                          for metric in distance_metrics}

# Function to calculate KDE and find the position of the max
def calculate_kde_max_position(data, bins):
    kde = gaussian_kde(data)
    kde_values = kde(bins)
    max_position = bins[np.argmax(kde_values)]
    return max_position

# Parameters for KDE and bootstrap
binsize = 0.25
bins = np.arange(0, 100, binsize)
bootstrap_samples = 5

# Initialize lists for storing results
distance_positions = []
distance_errors = []

# Bootstrap function
def bootstrap_error(data, csr, bins, n_samples):
    distances = []
    for _ in range(n_samples):
        resampled_data = resample(data)
        resampled_csr = resample(csr)
        max_data = calculate_kde_max_position(resampled_data, bins)
        max_csr = calculate_kde_max_position(resampled_csr, bins)
        
        # Calculate distance (CSR Max - Data Max)
        distance = max_csr - max_data
        distances.append(distance)
        
    return np.std(distances)

# Calculate distances and errors
for metric in distance_metrics:
    if len(combined_distances[metric]) > 0 and len(combined_distances_csr[metric]) > 0:
        # Calculate max positions
        max_position_data = calculate_kde_max_position(combined_distances[metric], bins)
        max_position_csr = calculate_kde_max_position(combined_distances_csr[metric], bins)
        
        # Calculate distance
        distance = max_position_csr - max_position_data
        
        # Calculate errors using bootstrapping
        dist_error = bootstrap_error(combined_distances[metric], combined_distances_csr[metric], bins, bootstrap_samples)
        
        # Store results
        distance_positions.append(distance)
        distance_errors.append(dist_error)

# Plot settings
colors = ['#2880C4', '#F4B942', '#D9534F', '#5CB85C']
labels = ['1NN', '2NN', '3NN', '4NN']

# Bar Plot for Distance with Error Bars and Caps
fig, ax = plt.subplots(figsize=(8, 6))

# Bar plot for distance with error bars (with caps)
bars = ax.bar(labels, distance_positions, yerr=distance_errors, color=colors, alpha=0.7, label='Distance', 
              capsize=5, error_kw=dict(lw=1.5))

# Set axis labels and title
ax.set_ylabel('Distance between Max Positions (CSR Max - Data Max) [nm]', color='black')
ax.tick_params(axis='y', labelcolor='black')
ax.set_title('Distance between Max Positions with Error Bars')

# Set y-axis limit slightly above the max error for better visibility
ax.set_ylim(0, max(np.array(distance_positions) + np.array(distance_errors)) * 1.2)

# Adjust layout and show plot
plt.tight_layout()
plt.show()

# Save the figure
fig.savefig(main_dir + '/distance_with_errors.pdf', format='pdf', bbox_inches='tight', pad_inches=0.1, transparent=True)





