#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 20:56:34 2023

@author: Luciano A. Masullo
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn

π = np.pi
plt.close('all')

experiments = ['MCF10A', 'MCF10AT', 'MCF10A (TGFβ)', 'MCF10AT (TGFβ)']

data = {}

data[experiments[0]] = np.array([13, 19, 18])
data[experiments[1]] = np.array([66, 45])
data[experiments[2]] = np.array([24, 32])    
data[experiments[3]] = np.array([8, 12, 10])


for i in range(4):
        
    plt.bar(experiments[i], np.mean(data[experiments[i]]), color='#2880C4',
            alpha=0.5, width=0.5, edgecolor='black', linewidth=1)
    
    # plt.errorbar(experiments[i], np.mean(data[experiments[i]]), 
    #              yerr=np.std(data[experiments[i]]), color='black', capsize=4)
    
    for j, value in enumerate(data[experiments[i]]):
        
        plt.scatter(experiments[i], value, color='gray', alpha=0.5)
        
plt.ylabel('Sugar density (μm^-2)')
               
