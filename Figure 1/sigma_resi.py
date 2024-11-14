#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 18:12:53 2024

@author: masullo
"""

import numpy as np
import matplotlib.pyplot as plt

# Data for the curve
K = np.linspace(1, 100, 500)  # Values of K from 1 to 100
sigma_localization = 3  # Arbitrary constant to control the curve scaling
sigma_resi = sigma_localization / np.sqrt(K)

# Plot
fig, ax = plt.subplots(figsize=(6, 5))
ax.plot(K, sigma_resi, color='navy', linewidth=3)

# Labels and styling
ax.set_xlabel(r'$K$', fontsize=16)
ax.set_ylabel(r'$\sigma_{\mathrm{RESI}} \, (\mathrm{nm})$', fontsize=16)
ax.set_xlim(1, 100)
ax.set_ylim(0, 3)

# Horizontal dotted line at sigma_resi = 1 nm
ax.axhline(y=1, color='black', linestyle='--', linewidth=1)
# Vertical dashed lines at K = 10 and K = 40
ax.axvline(x=10, color='black', linestyle='--', linewidth=0.8)
# ax.axvline(x=40, color='black', linestyle='--', linewidth=0.8)

# Annotation for the Å scale region
ax.fill_between(K, 0, 1, where=(K >= 10), color='gray', alpha=0.2)
ax.text(50, 0.5, r'Å scale', fontsize=14, ha='center', color='gray', fontweight='bold')

# Equation annotation
ax.text(30, 1.8, r'$\sigma_{\mathrm{RESI}} \approx \frac{\sigma_{\mathrm{Localization}}}{\sqrt{K}}$', 
        fontsize=18, ha='center')

# Show plot
plt.tight_layout()
plt.show()
