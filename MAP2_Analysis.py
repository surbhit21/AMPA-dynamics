#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:08:54 2021

@author: surbhitwagle
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py
import tables
import scipy.io

from PlotBinned import *
def MAP2_Analysis(data_file_path):
    molecule = "MAP2"
    dendritic_lens = np.array([25,50,75,100,125]); # dendritic lengths for which data is available
    analysis_for_lengths = np.array([0,0,0,1,0]); # dendritic lengths for which the analysis has to be done each value should be either 0 or 1
    count = np.array([0,0,0,0,0])
    step_size = 25;
    bin_size = 5; # for binned analysis
    scale = 1/4.8177; #um/pixel
    data_file = scipy.io.loadmat(data_file_path)
    data_dict = {};
    for lens in dendritic_lens:
        data_dict[int(lens/step_size)-1] = {}
    # data = np.array('data')
    # print(data_file['MAP2_data'].shape[1])
    for cells in data_file['MAP2_data'][0,3:]:
        for dendrite in cells[0,:]:
            dendrite_length = dendrite[0,2].shape[0]*scale;
            for lens in dendritic_lens:
                if dendrite_length > lens:
                    # print(dendrite[0,2].shape)
                    # print(int(lens/step_size))
                    cindex= int(lens/step_size) - 1
                    data_dict[int(lens/step_size)-1][count[cindex]] =  dendrite[0,2];
                    count[cindex]  += 1;
    # print(data_dict[0].keys())
    
    for i in range(0,len(analysis_for_lengths)):
        file_name = "MAP2_distribution_for_%i"%dendritic_lens[i];
        lab =  r'Mean $\pm$ S.E.M'
        y_lab = "Normalized Intensity"
        x_lab = r'Dendritic distance (in $\mu$M)'
        title = 'MAP2 intensity'
        if analysis_for_lengths[i] == 1:
            # print(count[i])
            print(file_name)
            PlotBinned(data_dict[i],dendritic_lens[i], scale, count[i], bin_size, file_name,title,lab,x_lab,y_lab, molecule,1)
if __name__ == '__main__':
    MAP2_Analysis("../../Neuron2019-data/Protein-data/MAP2_data.mat")