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

# Close all existing figures
plt.close('all')

# Define the main directories for ManNAz and GalNAz
dirs = {
    "ManNAz_combined": "/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz/ManNAz_combined",
    "GalNAz_combined": "/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/GalNAz/GalNAz_combined"
}

# Output directory
output_dir = "/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas"

# Initialize dictionaries to store the values
results = {
    "ManNAz_combined": {"clustered_fraction": [], "obs_density": []},
    "GalNAz_combined": {"clustered_fraction": [], "obs_density": []}
}

# Extract values from CSV files
for label, main_dir in dirs.items():
    for folder_name in sorted(os.listdir(main_dir)):
        folder_path = os.path.join(main_dir, folder_name)
        if os.path.isdir(folder_path) and folder_name.isdigit():
            results_path = os.path.join(folder_path, "integrated_results", "final_numbers.csv")
            if os.path.exists(results_path):
                df = pd.read_csv(results_path)
                results[label]["clustered_fraction"].append(df["clustered_fraction (rel. increase in %)"].iloc[0])
                results[label]["obs_density"].append(df["obs_density (μm^-2)"].iloc[0])

# Calculate mean and standard deviation
summary_stats = {}
for label in results.keys():
    cf_values = np.array(results[label]["clustered_fraction"])
    od_values = np.array(results[label]["obs_density"])
    summary_stats[label] = {
        "mean_cf": np.mean(cf_values),
        "std_cf": np.std(cf_values),
        "mean_od": np.mean(od_values),
        "std_od": np.std(od_values),
        "values_cf": cf_values,
        "values_od": od_values
    }

# Statistical tests
def significance_stars(p_value):
    return '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'n.s.'

# Two-sample t-tests between ManNAz and GalNAz
p_val_cf_groups = stats.ttest_ind(
    summary_stats["ManNAz_combined"]["values_cf"],
    summary_stats["GalNAz_combined"]["values_cf"]
).pvalue
p_val_od_groups = stats.ttest_ind(
    summary_stats["ManNAz_combined"]["values_od"],
    summary_stats["GalNAz_combined"]["values_od"]
).pvalue

# One-sample t-tests against zero
p_vals_zero = {
    "ManNAz_cf": stats.ttest_1samp(summary_stats["ManNAz_combined"]["values_cf"], 0).pvalue,
    "GalNAz_cf": stats.ttest_1samp(summary_stats["GalNAz_combined"]["values_cf"], 0).pvalue,
    "ManNAz_od": stats.ttest_1samp(summary_stats["ManNAz_combined"]["values_od"], 0).pvalue,
    "GalNAz_od": stats.ttest_1samp(summary_stats["GalNAz_combined"]["values_od"], 0).pvalue
}

# Significance stars
stars_cf_groups = significance_stars(p_val_cf_groups)
stars_od_groups = significance_stars(p_val_od_groups)

# Plotting
labels = ['ManNAz', 'GalNAz']
colors = ['#FF3C38', '#6C8EAD']
bar_width = 0.5

K = 40 # number of areas

# Clustered Fraction Plot
fig1, ax1 = plt.subplots(figsize=(4.5, 6))
means_cf = [summary_stats["ManNAz_combined"]["mean_cf"], summary_stats["GalNAz_combined"]["mean_cf"]]
stds_cf = [summary_stats["ManNAz_combined"]["std_cf"], summary_stats["GalNAz_combined"]["std_cf"]] / np.sqrt(K)

bars_cf = ax1.bar(labels, means_cf, yerr=stds_cf, color=colors, capsize=10, width=bar_width)
ax1.set_ylabel("Clustered Fraction (rel. increase in %)")
ax1.set_title("Comparison of Clustered Fraction")

# Annotate significance against zero
for i, (mean, std, key) in enumerate(zip(means_cf, stds_cf, ["ManNAz_cf", "GalNAz_cf"])):
    stars = significance_stars(p_vals_zero[key])
    ax1.text(i, mean + std * 1.1, stars, ha='center')

# Annotate comparison between groups
y_max_cf = max(means_cf) + max(stds_cf) * 1.4
ax1.plot([0, 1], [y_max_cf, y_max_cf], color='black')
ax1.text(0.5, y_max_cf * 1.02, stars_cf_groups, ha='center')

# Observed Density Plot
fig2, ax2 = plt.subplots(figsize=(4.5, 6))
means_od = [summary_stats["ManNAz_combined"]["mean_od"], summary_stats["GalNAz_combined"]["mean_od"]]
stds_od = [summary_stats["ManNAz_combined"]["std_od"], summary_stats["GalNAz_combined"]["std_od"]] / np.sqrt(K)

bars_od = ax2.bar(labels, means_od, yerr=stds_od, color=colors, capsize=10, width=bar_width)
ax2.set_ylabel("Observed Density (μm^-2)")
ax2.set_title("Comparison of Observed Density")

# Annotate significance against zero
for i, (mean, std, key) in enumerate(zip(means_od, stds_od, ["ManNAz_od", "GalNAz_od"])):
    stars = significance_stars(p_vals_zero[key])
    ax2.text(i, mean + std * 1.1, stars, ha='center')

# Annotate comparison between groups
y_max_od = max(means_od) + max(stds_od) * 1.4
ax2.plot([0, 1], [y_max_od, y_max_od], color='black')
ax2.text(0.5, y_max_od * 1.02, stars_od_groups, ha='center')

# Save figures as PDFs
fig1.savefig(os.path.join(output_dir, "Clustered_Fraction_Comparison.pdf"), bbox_inches='tight')
fig2.savefig(os.path.join(output_dir, "Observed_Density_Comparison.pdf"), bbox_inches='tight')

plt.show()






