#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:33:16 2022

@author: surbhitwagle
"""

win_len = 10
op_mat10 = GetSlidingWindowMeanMatrix(GFP_data[100],int(10/scale))
op_mat5 = GetSlidingWindowMeanMatrix(GFP_data[100],int(5/scale))
op_mat15 = GetSlidingWindowMeanMatrix(GFP_data[100],int(15/scale))
gfp_mean = GFP_data[100].mean(axis=0)
gfp_sw5_mean = op_mat5.mean(axis=0)
gfp_sw10_mean = op_mat10.mean(axis=0)
gfp_sw15_mean = op_mat15.mean(axis=0)
plt.plot(gfp_mean,label='m')
plt.plot(gfp_sw5_mean,label='5sw')
plt.plot(gfp_sw10_mean,label='10sw')
plt.plot(gfp_sw15_mean,label='15sw')
plt.legend()
plt.show()


win_len = 15
dend_len = 100
scale = 0.24
op_mat15 = GetSlidingWindowMeanMatrix(GFP_data[dend_len],int(win_len/scale))
op_surf_15 = GetSlidingWindowMeanMatrix(surf_glua2_data[dend_len],int(win_len/scale))
op_int_15 = GetSlidingWindowMeanMatrix(int_glua2_data[dend_len],int(win_len/scale))



gfp_sw15_mean = op_mat15.mean(axis=0)
surf_sw15_mean = op_surf_15.mean(axis=0)
int_sw15_mean = op_int_15.mean(axis=0)

surf_density = surf_sw15_mean/gfp_sw15_mean
int_density = int_sw15_mean/gfp_sw15_mean
surf_density = surf_density/surf_density[0]
int_density = int_density/int_density[0]
ratios = surf_density/int_density
plt.plot(bins[1:int(dend_len/scale)],surf_density[:int(dend_len/scale)-1],label='surf')
plt.plot(bins[1:int(dend_len/scale)],int_density[:int(dend_len/scale)-1],label='int')
plt.plot(bins[1:int(dend_len/scale)],ratios[:int(dend_len/scale)-1],label='ratio')
# plt.ylim([0.,2])
plt.legend()
plt.show()


