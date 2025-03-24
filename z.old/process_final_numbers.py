#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:18:35 2024

@author: masullo
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Specify the main directory path
# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Filapodia MCF10As/Filopodia'  # Adjust to your specific folder path if different
# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Filapodia MCF10As/Random areas' 

# main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Glycosylated spherical domains HMECs/Random FOVs'
main_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Glycosylated spherical domains HMECs/Spherical clusters'


# Initialize a list to store the extracted data
data = []

# Verify if the directory exists
if not os.path.exists(main_dir):
    print(f"{main_dir} does not exist.")
else:
    # Print the current working directory
    print("Current Working Directory:", os.getcwd())

    # Walk through the main directory
    for root, dirs, files in os.walk(main_dir):
        print("Current Directory:", root)  # Check what directories are being accessed
        
        # Strip whitespace from filenames in the list
        stripped_files = [file.strip() for file in files]
        
        # Check if the current path contains the final_numbers.csv file
        if "final_numbers.csv" in stripped_files:
            # Construct the full path to the file
            csv_path = os.path.join(root, "final_numbers.csv")
            print(f"Found CSV: {csv_path}")  # Confirm CSV file found

            try:
                # Read the file using pandas
                df = pd.read_csv(csv_path)
                print(f"CSV Headers for {csv_path}: {df.columns.tolist()}")  # Print headers
                
                # Check the shape of the DataFrame
                print(f"DataFrame shape for {csv_path}: {df.shape}")  # Print shape
                
                # Check if the DataFrame has at least one row
                if df.shape[0] > 0:  
                    # Extract the values from the first row (index 0)
                    values = df.iloc[0].tolist()
                    data.append(values)
                    print(f"Extracted values from {csv_path}: {values}")
                else:
                    print(f"Skipping {csv_path}: no data rows")
            
            except Exception as e:
                print(f"Error reading {csv_path}: {e}")

# Convert to a NumPy array
data_array = np.array(data)

# Calculate mean and standard deviation
means = np.mean(data_array, axis=0)
stds = np.std(data_array, axis=0)

# Round means and stds for densities to integers, others to 1 decimal place
means[0] = int(round(means[0]))  # Densities as integer
stds[0] = int(round(stds[0]))    # Densities as integer
means[1:] = np.round(means[1:], 1)  # Round other means to 1 decimal
stds[1:] = np.round(stds[1:], 1)    # Round other stds to 1 decimal

# Print results
print("Means:", means)
print("Standard Deviations:", stds)

# Save the data to a CSV file in the main_dir
output_csv_path = os.path.join(main_dir, "aggregated_data.csv")
output_df = pd.DataFrame(data_array, columns=df.columns.tolist())
output_df.to_csv(output_csv_path, index=False)
print(f"Saved aggregated data to {output_csv_path}")

# Plot histograms for each column
headers = df.columns.tolist()  # Get headers from the DataFrame for titles
plt.figure(figsize=(15, 5))

for i in range(data_array.shape[1]):
    plt.subplot(1, 3, i + 1)  # Create a subplot for each column
    plt.hist(data_array[:, i], bins=10, edgecolor='black', alpha=0.7)
    plt.title(headers[i])
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    
    # Format mean ± std based on the column
    mean_val = int(means[i]) if i == 0 else means[i]  # Convert to int if it's the density
    std_val = int(stds[i]) if i == 0 else stds[i]  # Round std for density
    
    if i == 0:
        # For density
        text = f'{mean_val} ± {std_val} μm⁻²'
    else:
        # For the other two columns
        text = f'{mean_val:.1f} ± {std_val:.1f} %'
    
    # Add formatted mean ± std to the top right of the histogram
    plt.text(0.65, 0.95, text, 
             transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

# Adjust y-limits for better visibility
for ax in plt.gcf().get_axes():
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)

plt.tight_layout()  # Adjust layout to prevent overlap

# Save the plot to a file in the main_dir as a PDF
plot_path = os.path.join(main_dir, "histograms.pdf")
plt.savefig(plot_path)
print(f"Saved histogram plot to {plot_path}")

plt.show()  # Display the histograms



