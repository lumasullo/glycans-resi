#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 18:05:42 2024

@author: masullo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

plt.close('all')

# Set parameters
mean1, mean2 = -0.4, 0.4       # Means for the two distributions (in nm)
sigma = 3                 # Standard deviation for both distributions (in nm)
N = 500                   # Number of samples per distribution

# Generate Gaussian-distributed samples
dist1 = np.random.normal(mean1, sigma, N)
dist2 = np.random.normal(mean2, sigma, N)

# Calculate the narrower sigma for the secondary Gaussian curves
sigma_narrow = sigma / np.sqrt(N)

# Define bin edges to ensure both histograms align
bin_edges = np.linspace(-10, 10, 30)

# Plot settings
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_box_aspect(1)  # Ensures a square aspect ratio for the plot area

# Colors and transparency
# color1, color2 = '#1D70B5', '#006633'  # Specified colors
color1, color2 = '#DF2935', '#3772FF'
alpha_val = 0.15  # Increased transparency level for histograms

# Plot histograms with Gaussian fits
ax.hist(dist1, bins=bin_edges, density=True, color=color1, alpha=alpha_val, edgecolor='black', linewidth=0.5, label='Mean 0 nm')
ax.hist(dist2, bins=bin_edges, density=True, color=color2, alpha=alpha_val, edgecolor='black', linewidth=0.5, label='Mean 1 nm')

# Gaussian fits
x_vals = np.linspace(-10, 10, 500)
fit1 = norm.pdf(x_vals, mean1, sigma)
fit2 = norm.pdf(x_vals, mean2, sigma)

# Plot Gaussian fits
ax.plot(x_vals, fit1, color=color1, linewidth=1.5, linestyle='--', label='Fit Mean 0 nm')
ax.plot(x_vals, fit2, color=color2, linewidth=1.5, linestyle='--', label='Fit Mean 1 nm')

# Narrower Gaussian curves with the same peak height as the broad Gaussians
narrow_fit1 = norm.pdf(x_vals, mean1, sigma_narrow)
narrow_fit2 = norm.pdf(x_vals, mean2, sigma_narrow)

# Scale the narrow Gaussians to match the peak height of the broad Gaussians
scaling_factor1 = fit1.max() / narrow_fit1.max()
scaling_factor2 = fit2.max() / narrow_fit2.max()

ax.plot(x_vals, narrow_fit1 * scaling_factor1, color=color1, linestyle='-', linewidth=2, label='Narrow Fit Mean 0 nm')
ax.plot(x_vals, narrow_fit2 * scaling_factor2, color=color2, linestyle='-', linewidth=2, label='Narrow Fit Mean 1 nm')

# Labels and legend
ax.set_title("Gaussian Distributions with Fits and Scaled Narrowed Curves")
ax.set_xlabel("Distance (nm)")
ax.set_ylabel("Density")

# ax.set_ylim(0.18)
# ax.legend()

plt.tight_layout()
plt.show()

