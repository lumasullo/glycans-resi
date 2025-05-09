#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 13:10:36 2025

@author: masullo
"""

import os
import shutil

# Define paths
base_folder = "/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/z.for_revisions/ManNAz/ManNAz_combined"
destination_folder = os.path.join(base_folder, "picks_combined")

# Create destination folder if it doesn't exist
os.makedirs(destination_folder, exist_ok=True)

# List of target files to process
target_files = ["target_picks.yaml", "target_picked.hdf5", "target_picked.yaml"]

# Loop through folders 1 to 40
for i in range(1, 41):
    folder_path = os.path.join(base_folder, str(i))
    
    for file_name in target_files:
        target_file = os.path.join(folder_path, file_name)

        # Check if the target file exists
        if os.path.exists(target_file):
            # Rename the file
            new_filename = f"{file_name.split('.')[0]}_{i}.{file_name.split('.')[-1]}"
            destination_path = os.path.join(destination_folder, new_filename)

            # Move and rename the file
            shutil.move(target_file, destination_path)
            print(f"Moved: {target_file} → {destination_path}")
        else:
            print(f"Warning: {target_file} not found!")

print("Task completed!")
