#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 18:46:48 2024

@author: masullo
"""

import numpy as np


def cart2sph(x, y, z):
    
   xy = np.sqrt(x**2 + y**2) # sqrt(x² + y²)
    
   x_2 = x**2
   y_2 = y**2
   z_2 = z**2

   r = np.sqrt(x_2 + y_2 + z_2) # r = sqrt(x² + y² + z²)

   theta = np.arctan2(y, x) 

   phi = np.arctan2(xy, z) 

   return r, theta, phi