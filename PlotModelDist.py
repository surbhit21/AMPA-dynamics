#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 18:16:38 2022

@author: surbhitwagle
"""

import numpy as np
ps_dist = np.load("Data/ps_dist.npy")
pc_dist = np.load("Data/pc_dist.npy")

from PlottingWidgetAMPA import *
pwa = PlottingWidgetAMPA()
x = np.arange(0,500,0.24)
pwa.PlotSingleSimTwoProtein(x, ps_dist, pc_dist, "surface", "cytoplasmic", r"Dendritic distance (in \mu M)", "GluA2 copy number", "Steady state distribution", "fitted_dist_GluA2")