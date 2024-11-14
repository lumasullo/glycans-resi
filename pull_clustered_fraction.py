#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 17:33:02 2024

@author: masullo
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# Define the main directories for ManNAz and GalNAz
dirs = {
    "ManNAz_combined": "/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz/ManNAz_combined",
    "GalNAz_combined": "/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/GalNAz/GalNAz_combined"
}

# Initialize dictionaries to store the values
results = {
    "ManNAz_combined": {"clustered_fraction": [], "obs_density": []},
    "GalNAz_combined": {"clustered_fraction": [], "obs_density": []}
}

# Iterate through each directory and extract the values
for label, main_dir in dirs.items():
    for folder_name in sorted(os.listdir(main_dir)):
        folder_path = os.path.join(main_dir, folder_name)

        # Check if the folder is a directory and named with a number
        if os.path.isdir(folder_path) and folder_name.isdigit():
            results_path = os.path.join(folder_path, "integrated_results", "final_numbers.csv")

            if os.path.exists(results_path):
                df = pd.read_csv(results_path)

                if "clustered_fraction (rel. increase in %)" in df.columns:
                    clustered_fraction_value = df["clustered_fraction (rel. increase in %)"].iloc[0]
                    results[label]["clustered_fraction"].append(clustered_fraction_value)

                if "obs_density (μm^-2)" in df.columns:
                    obs_density_value = df["obs_density (μm^-2)"].iloc[0]
                    results[label]["obs_density"].append(obs_density_value)

# Calculate the mean and standard deviation for both metrics
summary_stats = {}
for label in results.keys():
    clustered_fraction_values = np.array(results[label]["clustered_fraction"])
    obs_density_values = np.array(results[label]["obs_density"])

    summary_stats[label] = {
        "mean_clustered_fraction": np.mean(clustered_fraction_values),
        "std_clustered_fraction": np.std(clustered_fraction_values),
        "mean_obs_density": np.mean(obs_density_values),
        "std_obs_density": np.std(obs_density_values),
        "values_clustered_fraction": clustered_fraction_values,
        "values_obs_density": obs_density_values
    }

# Statistical testing (t-test)
def significance_stars(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'n.s.'

# T-test for clustered_fraction
t_stat_cf, p_val_cf = stats.ttest_ind(
    summary_stats["ManNAz_combined"]["values_clustered_fraction"],
    summary_stats["GalNAz_combined"]["values_clustered_fraction"]
)
stars_cf = significance_stars(p_val_cf)

# T-test for obs_density
t_stat_od, p_val_od = stats.ttest_ind(
    summary_stats["ManNAz_combined"]["values_obs_density"],
    summary_stats["GalNAz_combined"]["values_obs_density"]
)
stars_od = significance_stars(p_val_od)

# Plotting
labels = ['ManNAz', 'GalNAz']
colors = ['#FF3C38', '#6C8EAD']  # Red for ManNAz, Blue-gray for GalNAz
bar_width = 0.5

# Figure 1: Clustered Fraction Comparison
fig1, ax1 = plt.subplots(figsize=(7, 7))
means_clustered_fraction = [summary_stats["ManNAz_combined"]["mean_clustered_fraction"],
                            summary_stats["GalNAz_combined"]["mean_clustered_fraction"]]
stds_clustered_fraction = [summary_stats["ManNAz_combined"]["std_clustered_fraction"],
                           summary_stats["GalNAz_combined"]["std_clustered_fraction"]]

ax1.bar(labels, means_clustered_fraction, yerr=stds_clustered_fraction, color=colors, capsize=10, width=bar_width)
ax1.set_ylabel("Clustered Fraction (rel. increase in %)")
# ax1.set_title("Comparison of Clustered Fraction")
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Adjust y-axis limits
ax1.set_ylim(0, max(means_clustered_fraction) + max(stds_clustered_fraction) * 2)

# Add significance stars
y_max_cf = max(means_clustered_fraction) + max(stds_clustered_fraction)
ax1.text(0.5, y_max_cf * 1.1, stars_cf, fontsize=15, ha='center', color='black')
ax1.plot([0, 1], [y_max_cf * 1.05, y_max_cf * 1.05], color='black', linewidth=1.5)

# Figure 2: Observed Density Comparison
fig2, ax2 = plt.subplots(figsize=(7, 7))
means_obs_density = [summary_stats["ManNAz_combined"]["mean_obs_density"],
                     summary_stats["GalNAz_combined"]["mean_obs_density"]]
stds_obs_density = [summary_stats["ManNAz_combined"]["std_obs_density"],
                    summary_stats["GalNAz_combined"]["std_obs_density"]]

ax2.bar(labels, means_obs_density, yerr=stds_obs_density, color=colors, capsize=10, width=bar_width)
ax2.set_ylabel("Observed Density (μm^-2)")
# ax2.set_title("Comparison of Observed Density")
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Adjust y-axis limits
ax2.set_ylim(0, max(means_obs_density) + max(stds_obs_density) * 2)

# Add significance stars
y_max_od = max(means_obs_density) + max(stds_obs_density)
ax2.text(0.5, y_max_od * 1.1, stars_od, fontsize=15, ha='center', color='black')
ax2.plot([0, 1], [y_max_od * 1.05, y_max_od * 1.05], color='black', linewidth=1.5)

# Show the plots
plt.tight_layout()
plt.show()



