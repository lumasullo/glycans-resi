#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:20:23 2024

@author: masullo
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# Set the main directory
main_dir = '/Volumes/pool-miblab/users/masullo/z.fromKareem/main experiments/Glycosylated spherical domains HMECs'
# main_dir = '/Volumes/pool-miblab/users/masullo/z.fromKareem/main experiments/Filapodia MCF10As'

# Dynamically find folders inside the main directory
directories = [os.path.join(main_dir, folder) for folder in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, folder))]

# Check if there are exactly two folders
if len(directories) != 2:
    raise ValueError(f"Expected exactly 2 folders, but found {len(directories)}: {directories}")

# Initialize a dictionary to store data for each category
data_dict = {'Density': [], 'Clustered Fraction': [], 'NND Peak Shift': []}

# Loop through the directories and load the aggregated_data.csv files
for dir_path in directories:
    csv_path = os.path.join(dir_path, 'aggregated_data.csv')
    
    if os.path.exists(csv_path):
        # Load the CSV file into a DataFrame
        df = pd.read_csv(csv_path)

        # Strip whitespace from column names
        df.columns = df.columns.str.strip()

        # Print the headers for debugging
        print(f"Headers in {csv_path}: {df.columns.tolist()}")  # Debugging line
        
        # Check if the expected columns exist before appending
        if 'obs_density (μm^-2)' in df.columns:
            data_dict['Density'].append(df['obs_density (μm^-2)'].values)
        else:
            print(f"'obs_density (μm^-2)' not found in {csv_path}")
        
        if 'clustered_fraction (rel. increase in %)' in df.columns:
            data_dict['Clustered Fraction'].append(df['clustered_fraction (rel. increase in %)'].values)
        else:
            print(f"'clustered_fraction (rel. increase in %)' not found in {csv_path}")
        
        if 'nnd_peak_shift (nm)' in df.columns:
            data_dict['NND Peak Shift'].append(df['nnd_peak_shift (nm)'].values)
        else:
            print(f"'nnd_peak_shift (nm)' not found in {csv_path}")
    else:
        print(f"{csv_path} does not exist.")

# Convert lists to arrays for boxplotting and check for empty data
for key in data_dict.keys():
    data_dict[key] = np.array(data_dict[key], dtype=object)  # Convert to a 2D array

# Create box plots for each category
fig, axs = plt.subplots(1, 3, figsize=(12, 5))  # Adjust the figure size

light_gray = '#D3D3D3'  # Light gray color for box plots

# Define y-axis labels for each subplot
y_labels = ['Density (μm⁻²)', 'Rel. change (%)', 'Rel. change (%)']

for i, key in enumerate(data_dict.keys()):
    # Ensure there is data to plot
    if len(data_dict[key]) > 0 and all(len(data) > 0 for data in data_dict[key]):
        # Create a box plot for each category
        box = axs[i].boxplot(data_dict[key], labels=[os.path.basename(directories[j]) for j in range(len(data_dict[key]))], patch_artist=True)

        # Fill boxes with light gray color
        for patch in box['boxes']:
            patch.set_facecolor(light_gray)
            patch.set_edgecolor('black')
            patch.set_linewidth(1)

        axs[i].set_title(key)
        axs[i].set_ylabel(y_labels[i])  # Set the y-axis label
        axs[i].grid(axis='y', linestyle='--', alpha=0.7)

        # Set the median line color to black
        for median in box['medians']:
            median.set_color('black')

        # Perform ANOVA to test for significance
        f_val, p_val = stats.f_oneway(*data_dict[key])

        # Prepare significance stars based on p-value
        def significance_stars(p_value):
            if p_value < 0.001:
                return '***'
            elif p_value < 0.01:
                return '**'
            elif p_value < 0.05:
                return '*'
            else:
                return 'ns'  # not significant

        # Get stars for the plot
        stars = significance_stars(p_val)

        # Annotate significance on the plot
        y_max = np.max(np.concatenate(data_dict[key])) * 0.95  # Get max y value for annotation
        axs[i].text(1.5, y_max, stars, fontsize=12, ha='center', color='black')

        # Draw a line for significance
        axs[i].plot([1, 2], [y_max, y_max], color='black', linestyle='-', linewidth=1.5)

    else:
        axs[i].set_title(f"{key} - No Data")
        axs[i].axis('off')  # Turn off axis if there's no data

# Adjust layout for tighter spacing
plt.subplots_adjust(wspace=0.1)  # Adjust the width space between subplots
plt.tight_layout(pad=2.0)  # Adjust padding for a tighter layout
plt.show()



