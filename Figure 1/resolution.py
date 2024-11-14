#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 18:51:41 2024

@author: masullo
"""

import matplotlib.pyplot as plt
import numpy as np

# Data for microscopy methods and their resolutions
methods = ['Diff. limited', 'STORM', 'DNA-PAINT', 'RESI']
resolutions = [260, 25, 8, 0.9]

# Colors for each bar (darker for lower resolutions)
colors = ['#b3cde0', '#6497b1', '#005b96', '#03396c']

# Plotting the horizontal bar chart
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.barh(methods, resolutions, color=colors, edgecolor='black', height=0.6)

# Set the x-axis to a log scale for better differentiation of smaller values
ax.set_xscale('log')

# ax.set_aspect(1)

# Add annotations on the bars
for bar, res in zip(bars, resolutions):
    ax.text(
        bar.get_width() + 2,  # Position right of the bar
        bar.get_y() + bar.get_height() / 2,
        f"{res} nm",
        va='center',
        ha='left',
        color='black',
        fontsize=12
    )

# Title and labels
ax.set_title('Resolution Comparison of Microscopy Methods', fontsize=16, fontweight='bold')
ax.set_xlabel('Resolution (nm)', fontsize=14)
ax.invert_yaxis()  # Invert y-axis for ranking order

# Customizing the x-axis for better readability
ax.grid(visible=True, which="both", linestyle='--', color='gray', alpha=0.5)
ax.set_axisbelow(True)

plt.tight_layout()
plt.show()
