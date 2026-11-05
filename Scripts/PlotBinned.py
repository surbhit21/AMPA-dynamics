#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:09:23 2021

@author: surbhitwagle
"""

import matplotlib.pyplot as plt
import numpy as np
from math import ceil
from numpy.polynomial.polynomial import polyfit
import os
from pathlib import Path
def PlotBinned(data_cells,dendritic_len,scale,count,bin_size,file_name,title_string,label_string,x_lab,y_lab,molecule,save_it,im_width=10,im_height=6):
    fig, ax = plt.subplots(figsize=(im_width, im_height))
    plt.rc('font', **{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    # print(data_cells)
    bins = np.arange(0,dendritic_len,bin_size);
    binned_sum = np.zeros((count,len(bins)))
    norm_binned_sum = np.zeros((count,len(bins)))
    for i in range(0,count):
        x = data_cells[i][0:ceil(dendritic_len/scale),0]
        data = data_cells[i][0:ceil(dendritic_len/scale),1];
        # print(data)
        digitized = np.digitize(x, bins)
        # print(digitized)
        data_bin_sum = [data[digitized == j].sum() for j in range(1, len(bins)+1)]
        binned_sum[i] = data_bin_sum
        norm_binned_sum[i] = data_bin_sum/data_bin_sum[0]
    mean_norm_binned = norm_binned_sum.mean(0)
    std_norm_binned = norm_binned_sum.std(0)
    sem_norm_binned = std_norm_binned/np.sqrt(count)
    plt.errorbar(bins+bin_size/2,mean_norm_binned,sem_norm_binned,label=label_string,marker='d')
    plt.title(title_string)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    # Fit with polyfit
    x = np.arange(0,dendritic_len,scale)
    b, m = polyfit(bins+bin_size/2, mean_norm_binned, 1)
    plt.plot(x, b + m * x, '-',label='linear-fit')
    per_drop_at_len = -1*m*100*100;
    title_string += '\n'+r' bin size = %i $\mu$M, %s drop in intesnity at lenght %i $\mu$M is %.2f'%(bin_size,'%',dendritic_len,per_drop_at_len)
    plt.legend()
    folder = os.path.dirname(os.path.abspath(__file__)) + "/Figures/%s"%(molecule)
    Path(folder).mkdir(parents=True, exist_ok=True)
    if save_it == 1:
        plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=150)
        plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=150)
        print("saved figures to: {%s/%s}" %(folder, file_name))
    else:
        print("Plots not saved")
    plt.show()
    # print(mean_norm_binned,sem_norm_binned)