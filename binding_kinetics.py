#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 12:41:33 2025

@author: masullo
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import groupby
from operator import itemgetter
import os
import yaml

def load_info(path):
    path_base, path_extension = os.path.splitext(path)
    filename = path_base + ".yaml"
    try:
        with open(filename, "r") as info_file:
            info = list(yaml.load_all(info_file, Loader=yaml.FullLoader))
    except FileNotFoundError:
        print("\nAn error occured. Could not find metadata file:\n{}".format(filename))
    return info

def BindingEvents(frames):
    
    ranges =[]
    
    for k, g in groupby(enumerate(frames),lambda x:x[0]-x[1]):
                
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        ranges.append((group[0],group[-1]))
        
    events = np.asarray(ranges)
    beginnings = (events[:,0])
    dark_times = beginnings[1:] - beginnings[:-1]
    bright_times = events[:,1] - events[:,0] + 1
    
    return (events, bright_times, dark_times)


plt.close('all')

# Define the path and filenames
path = '/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/z.tests_binding_kinetics/HMEC 240617/'
filenames = [
    'R1_125pM_35mW_1_MMStack_Pos0.ome_locs_render_aligned_apicked_filter_clustered.hdf5',
    'R2_50pM_35mW_1_MMStack_Pos0.ome_locs_render_aligned_apicked_filter_clustered.hdf5',
    'R3_50pM_35mW_1_MMStack_Pos0.ome_locs_render_aligned_apicked_filter_clustered.hdf5',
    'R4_5pM_35mW_1_MMStack_Pos0.ome_locs_render_aligned_apicked_filter_clustered.hdf5',
    'R5_125pM_35mW_1_MMStack_Pos0.ome_locs_render_aligned_apicked_filter_clustered.hdf5',
    'R6_125pM_35mW_1_MMStack_Pos0.ome_locs_render_aligned_apicked_filter_clustered.hdf5'
]

exposure_time = 0.1 # in seconds

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, filename in enumerate(filenames):
    filepath = path + filename
    
    info = load_info(filepath)
    nframes = info[0]['Frames']
    df = pd.read_hdf(filepath, key='locs')
    
    cluster_iter = np.sort(np.unique(np.array(df['group'])))
    n_binding_events = []
    dark_times = []
    bright_times = []
    
    for i in cluster_iter:
        
        print(i)
        
        cluster = df[df['group'] == i]
        on_frames = np.array(cluster['frame'])
        events, bt, dt = BindingEvents(on_frames)
        n_binding_events.append(len(events))
        dark_times.append(np.mean(dt))
        bright_times.append(np.mean(bt))
    
    n_binding_events = np.array(n_binding_events)
    n_binding_events = n_binding_events[n_binding_events < 50] # filter for outliers
    
    axes[idx].hist(n_binding_events, bins=range(0, 50), alpha=0.5, density=True, 
                   color='#3767E1', edgecolor='black', linewidth=0.5)
    axes[idx].set_title(f'R{idx + 1}')
    axes[idx].set_xlabel('Number of Binding Events')
    axes[idx].set_ylabel('Density')
    
    axes[idx].set_xlim(0, 40)
    axes[idx].set_ylim(0, 0.23)
    
    axes[idx].tick_params(direction='in', length=2, width=1)
    
    n_mean, n_std = np.mean(n_binding_events), np.std(n_binding_events)    
    n_median = np.median(n_binding_events)
    
    n_mean = np.around(n_mean, 1)
    n_median = np.around(n_median, 1)
    n_std = np.around(n_std, 1)
    
    axes[idx].text(20, 0.12, f'n = {n_mean:.1f} ± {n_std:.1f}')
    
    bright_times = np.array(bright_times)
    bt_mean = bright_times.mean()
    bt_std = bright_times.std()
    
    bt_mean = np.around(bt_mean, 2) * exposure_time
    bt_std = np.around(bt_std, 2) * exposure_time
    
    axes[idx].text(20, 0.10, f'bt = ({bt_mean:.2f} ± {bt_std:.2f}) s')
    
    locs_mean = n_mean * (bt_mean / exposure_time)
    locs_std = n_std * (bt_mean / exposure_time)
    
    axes[idx].text(20, 0.08, f'n_locs = {locs_mean:.1f} ± {locs_std:.1f}')

    
    
plt.tight_layout()
plt.show()

print('')

# -------- Save Source Data --------

# Collect per-replicate results
replicate_data = []

for idx, filename in enumerate(filenames):
    replicate = f'R{idx + 1}'
    
    filepath = path + filename
    info = load_info(filepath)
    df = pd.read_hdf(filepath, key='locs')
    
    cluster_iter = np.sort(np.unique(np.array(df['group'])))
    per_cluster_events = []
    
    for i in cluster_iter:
        cluster = df[df['group'] == i]
        on_frames = np.array(cluster['frame'])
        events, bt, dt = BindingEvents(on_frames)
        if len(events) < 50:
            replicate_data.append({
                'replicate': replicate,
                'filename': filename,
                'cluster_id': int(i),
                'n_binding_events': len(events),
                'bright_time_mean': np.mean(bt) * exposure_time,
                'dark_time_mean': np.mean(dt) * exposure_time
            })

# Save to CSV
source_df = pd.DataFrame(replicate_data)
source_csv_path = os.path.join(path, 'EDF6_SourceData.csv')
source_df.to_csv(source_csv_path, index=False)
print(f'\nSource data saved to: {source_csv_path}')


# -------- Replot from Saved Source Data --------

df_loaded = pd.read_csv(source_csv_path)

fig, axes2 = plt.subplots(2, 3, figsize=(15, 10))
axes2 = axes2.flatten()

for idx in range(6):
    rep = f'R{idx + 1}'
    sub_df = df_loaded[df_loaded['replicate'] == rep]
    
    n_events = sub_df['n_binding_events']
    bt_mean = sub_df['bright_time_mean'].mean()
    bt_std = sub_df['bright_time_mean'].std()
    
    n_mean = n_events.mean()
    n_std = n_events.std()
    
    ax = axes2[idx]
    ax.hist(n_events, bins=range(0, 50), alpha=0.5, density=True,
            color='#3767E1', edgecolor='black', linewidth=0.5)
    
    ax.set_title(rep)
    ax.set_xlabel('Number of Binding Events')
    ax.set_ylabel('Density')
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 0.23)
    ax.tick_params(direction='in', length=2, width=1)
    
    ax.text(20, 0.12, f'n = {n_mean:.1f} ± {n_std:.1f}')
    ax.text(20, 0.10, f'bt = ({bt_mean:.2f} ± {bt_std:.2f}) s')
    
    locs_mean = n_mean * (bt_mean / exposure_time)
    locs_std = n_std * (bt_mean / exposure_time)
    ax.text(20, 0.08, f'n_locs = {locs_mean:.1f} ± {locs_std:.1f}')

plt.tight_layout()
plt.show()