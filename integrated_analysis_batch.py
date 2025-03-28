#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 12:40:39 2024

@author: masullo
"""

import os
import pandas as pd

# Define the root directory where the folders are stored

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/ManNAz'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/GalNAz'

# root_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/HMECs Homogenous areas/ManNAz'

# root_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/z.for_revisions/ManNAz'

root_dir = '/Users/masullo/Library/CloudStorage/Dropbox/z.forKareem_datashare/07.data_sharing/2024/Paper/z.for_revisions/GalNAz'



# Load the external Python script as a string
analysis_script_path = '/Users/masullo/Documents/GitHub/glycans-resi/integrated_analysis.py'
with open(analysis_script_path, 'r') as file:
    analysis_script = file.read()

# Loop through the main folders like "Filopodia", "Protrusions", etc.
for main_folder in os.listdir(root_dir):
    main_folder_path = os.path.join(root_dir, main_folder)

    # Check if the current item is a folder
    if os.path.isdir(main_folder_path):
        print(f'Processing folder: {main_folder}')
        
        # Loop through the subfolders (e.g., 1, 2, 3, ...)
        for subfolder in os.listdir(main_folder_path):
            subfolder_path = os.path.join(main_folder_path, subfolder)
            
            # Skip .DS_Store and non-directory items
            if subfolder == '.DS_Store' or not os.path.isdir(subfolder_path):
                continue
            
            # Construct the file path to the "target_picked.hdf5" file
            target_file_path = os.path.join(subfolder_path, 'target_picked.hdf5')

            # Check if the target file exists
            if os.path.isfile(target_file_path):
                print(f'Processing file: {target_file_path}')
                
                # Load the HDF5 file into a pandas DataFrame
                try:
                    df_exp = pd.read_hdf(target_file_path, key='locs')

                    # Inject the loaded DataFrame and file path into the script
                    exec(analysis_script, globals(), {'df_exp': df_exp, 'filepath': target_file_path})
                    
                    print(f'Successfully ran analysis for {target_file_path}')
                except Exception as e:
                    print(f'Error processing {target_file_path}: {e}')
            else:
                print(f'File not found: {target_file_path}')