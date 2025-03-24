#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 13:03:24 2025

@author: masullo
"""

import numpy as np
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import configparser
from timeit import default_timer as timer
from scipy.stats import gaussian_kde

from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import pdist


abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

from Functions import io
from Functions import dbscan

plt.close('all')

# # 240617_HMEC ManNAz
# path = r'/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/240617_HMEC/RESI_ManNAz/workflow_analysis/00_cluster_241009-1223/00_cluster_aggregation_241009-1223/04_save_datasets_aggregated/integrated_results/'
# filename = r'target_picked_dbscan_10_1.hdf5'

# 240618_HMEC GalNAz
path = r'/Volumes/pool-miblab/users/masullo/z_raw/GlycoRESI/240618_HMEC/RESI_GalNAz/workflow_analysis/00_cluster_241007-1052/00_cluster_aggregation_241007-1052/04_save_datasets_aggregated/integrated_results/'
filename = r'target_picked_dbscan_10_1.hdf5'

filepath = path + filename # comment / uncomment this line for batch analysis or single analysis respectively

df = pd.read_hdf(filepath, key = 'locs')

df_0 = df[df['nsugars_in_cluster'] == 1]
df_1 = df[df['nsugars_in_cluster'] > 1]

info = io.load_info(filepath.replace('.hdf5', f'.yaml'))

df_0_ = df_0.to_records() # convert to recarray

# save locs in Picasso file
df_0_filename = '%s_nonclustered.hdf5' % (filename.replace('.hdf5', ''))
io.save_locs(path + df_0_filename, df_0_, info)

df_1_ = df_1.to_records() # convert back to recarray

# save locs in Picasso file
df_1_filename = '%s_clustered.hdf5' % (filename.replace('.hdf5', ''))
io.save_locs(path + df_1_filename, df_1_, info)



