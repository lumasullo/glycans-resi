#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 12:40:39 2024

@author: masullo
"""

import os
import pandas as pd

# Define the root directory where the folders (e.g., "Filopodia" and "Protrusions") are stored

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/Advanced ROI to analyze/240919_MCF10A'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Full cells analysis/MCF10As'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Full cells analysis/HMECS'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Full cells analysis/aux_remaining_0'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Full cells analysis/aux_remaining_1'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Full cells analysis/aux_remaining_2'



# root_dir = '//Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Filapodia MCF10As/Random areas'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Filapodia MCF10As/Filopodia'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Glycosylated spherical domains HMECs/Random FOVs'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/main experiments/Glycosylated spherical domains HMECs/Spherical clusters'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/spherical nanodomains'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/Homogenous areas/ManNAz'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/High density nanoclusters/GalNAz'
# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/High density nanoclusters/ManNAz'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/MCF10A_and_MCF10AT/Spherical MCF10As'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/MCF10A_and_MCF10AT/Homogenous  ROIs'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/diff_density_areas/MCF10A_and_MCF10AT/High density domains'

# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/ManNAz'

root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/HMECs Homogenous areas/GalNAz'


# root_dir = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.fromKareem/homogenous areas data/MCF10As Homogenous  ROIs/'

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