# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 11:12:22 2022

@author: reinhardt
"""
import yaml
import os
import h5py
import numpy as np


# adapted from Picasso io.py
def load_info(path):
    path_base, path_extension = os.path.splitext(path)
    filename = path_base + ".yaml"
    try:
        with open(filename, "r") as info_file:
            info = list(yaml.load_all(info_file, Loader=yaml.FullLoader))
    except FileNotFoundError:
        print("\nAn error occured. Could not find metadata file:\n{}".format(filename))
    return info


# adapted from Picasso io.py
def save_info(path, info, default_flow_style=False):
    with open(path, "w") as file:
        yaml.dump_all(info, file, default_flow_style=default_flow_style)
        
    
# adapted from Picasso io.py
def save_locs(path, locs, info):
    #locs = _lib.ensure_sanity(locs, info)
    with h5py.File(path, "w") as locs_file:
        locs_file.create_dataset("locs", data=locs)
    base, ext = os.path.splitext(path)
    info_path = base + ".yaml"
    save_info(info_path, info)
    
def save_pos(output_path, filename, width, height, pos, info):
    # Save coordinates to csv or Picasso hdf5 file

    # filename without '.hdf5'
    frames = np.full(len(pos), int(0))
    x = pos[:, 0]
    y = pos[:, 1]
    lpx = np.full(len(pos), 0.001) # Dummy value required for Picasso Render to display points
    lpy = np.full(len(pos), 0.001)     
    photons = np.full(len(pos), 1)
    sx = np.full(len(pos), 0.001)  
    sy = np.full(len(pos), 0.001)  
    bg = np.full(len(pos), 0)

    
    LOCS_DTYPE = [
         ("frame", "u4"),
         ("x", "f4"),
         ("y", "f4"),
         ("lpx", "f4"),
         ("lpy", "f4"),
         ("photons", "f4"),
         ("sx", "f4"),
         ("sy", "f4"),
         ("bg", "f4")
         ]
    
    locs = np.rec.array(
         (frames, x, y, lpx, lpy, photons, sx, sy, bg),
         dtype=LOCS_DTYPE,
         )

    save_locs(os.path.join(output_path, filename + '.hdf5'), locs, info)
    
# Function to extract the 'Total Picked Area (um^2)' from the yaml file
def extract_total_picked_area(file_path):
    with open(file_path, 'r') as file:
        documents = yaml.safe_load_all(file)  # Use safe_load_all to handle multiple documents

        # Iterate over all documents in the YAML file
        for doc in documents:
            # Extract the 'Total Picked Area (um^2)' if it's found in the current document
            if 'Total Picked Area (um^2)' in doc:
                return doc['Total Picked Area (um^2)']

    return None  # Return None if the parameter is not found