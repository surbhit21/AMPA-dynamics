#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:26:13 2022

@author: surbhitwagle
"""
import math
import os.path
import argparse
from AMPA_model import RunModelWithFile
from GetPublishedData import *

import math
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
import scipy.stats as sst
from TemporalIntegration import baseline_param_file, dt
from Utility import *

"""
Plotting the complete profile dynamics over simulated time
"""
# fig,ax = plt.subplots()
# title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# line, = ax.plot(x_grid,  data[labels[0]][:,0], color = COLORS_dict["shaft_i"], lw=1,label=labels[0])
# line2, = ax.plot(x_grid, data[labels[1]][:,0], color = COLORS_dict["shaft_s"], lw=1,label=labels[1])
# line3, = ax.plot(x_grid, data[labels[2]][:,0], color = COLORS_dict["spine_s"], lw=1,label=labels[2])
#
# # txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# plt.legend(loc="upper right")
# plt.ylim([0,60])
# def animate(i):
#     line.set_data(x_grid,  data[labels[0]][:,i])
#     line2.set_data(x_grid,  data[labels[1]][:,i])
#     line3.set_data(x_grid,  data[labels[2]][:,i])
#     title.set_text('time = %.3f mins'  % (timepoints[i]/60))
#     # txt1.set_text(r"stim location = {} $\mu m$".format(locn))
#     return  line, line2, line3, title
#
# ani = FuncAnimation(fig, animate, interval=interval, blit=True, repeat=True, frames=f_rames)
# ani.save("{}/{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()

"""
Plotting a time gif of ratio of change from initial state
"""

# fig,ax = plt.subplots()
# line, = ax.plot(x_grid,  data[labels[0]][:,0]/data[labels[0]][:,0], color = COLORS_dict["shaft_i"], lw=1,label=labels[0])
# line2, = ax.plot(x_grid, data[labels[1]][:,0]/data[labels[1]][:,0], color = COLORS_dict["shaft_s"], lw=1,label=labels[1])
# line3, = ax.plot(x_grid, data[labels[2]][:,0]/data[labels[2]][:,0], color = COLORS_dict["spine_s"], lw=1,label=labels[2])
# # txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# plt.legend(loc="upper left")
# plt.plot(x_grid,  np.ones(x_grid.shape), color = 'black', lw=1,label=labels[0])
# plt.ylim([0,1.3])
# title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# def animateRatio(i):
#     line.set_data(x_grid,  data[labels[0]][:,i]/data[labels[0]][:,0])
#     line2.set_data(x_grid,  data[labels[1]][:,i]/data[labels[1]][:,0])
#     line3.set_data(x_grid,  data[labels[2]][:,i]/data[labels[2]][:,0])
#     title.set_text('time = %.3f mins'  % (timepoints[i]/60))
#     # txt1.set_text(r"stim location = {} $\mu m$".format(locn))
#     return line, line2, line3, title
#
# aniratio = FuncAnimation(fig, animateRatio, interval=interval, blit=True, repeat=True, frames=f_rames)
# aniratio.save("{}/ratio_{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()


"""
Plotting a time gif of zoomed at location of stimlation loc with a surrounding area of span 
"""
# fig,ax = plt.subplots()
# # loc =  50
# span = 50
# index_min = int((locn - span)/dx)
# index_max = int((locn + span)/dx)
# # beta_minx = int((locn-beta_span)/dx)
# # beta_maxx = int((locn+beta_span)/dx)
# # alpha_minx = int((locn-alpha_span)/dx)
# # alpha_maxx = int((locn+alpha_span)/dx)
# # range_new = [index_min:index_max]
# x_grid_new = x_grid[index_min:index_max]
# title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# line, = ax.plot(x_grid_new,  data[labels[0]][index_min:index_max,0], color = 'blue', lw=1,label=labels[0])
# line2, = ax.plot(x_grid_new, data[labels[1]][index_min:index_max,0], color = 'red', lw=1,label=labels[1])
# line3, = ax.plot(x_grid_new, data[labels[2]][index_min:index_max,0], color = 'green', lw=1,label=labels[2])
# ax.plot(x_grid_new,np.ones(x_grid_new.shape),"k--")
# txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# # min_line = ax.vlines(x= x_grid[beta_minx],ymin=0.9,ymax=1.1,linestyles="dashed")
# # max_line = ax.vlines(x= x_grid[beta_maxx],ymin=0.9,ymax=1.1,linestyles="dashed")
# # min_line = ax.vlines(x= x_grid[alpha_minx],ymin=0.9,ymax=1.1,linestyles="-",color='k')
# # max_line = ax.vlines(x= x_grid[alpha_maxx],ymin=0.9,ymax=1.1,linestyles="-",color="k")
# plt.legend(loc="upper right")
# plt.ylim([0.5,1.4])
#
# def animateZoomed(i):
#     line.set_data(x_grid_new, data[labels[0]][index_min:index_max, i]/data[labels[0]][index_min:index_max, 0])
#     line2.set_data(x_grid_new, data[labels[1]][index_min:index_max, i]/data[labels[1]][index_min:index_max, 0])
#     line3.set_data(x_grid_new, data[labels[2]][index_min:index_max, i]/data[labels[2]][index_min:index_max, 0])
#     title.set_text('time = %.3f mins' % ((timepoints[i]-30)/60))
#     txt1.set_text(r"stim location = {} $\mu m$".format(locn))
#     return line, line2, line3, title, txt1
#
# anizoomed = FuncAnimation(fig, animateZoomed, interval=interval, blit=True, repeat=True, frames=f_rames)
# anizoomed.save("{}/zoomed_{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()

"""
Plotting the time evolution of stimulated locations over time
"""
# print(timepoints.shape)
# breakpoint()
patterson_1g = patterson_data_1g()
patterson_4c = patterson_data_4c()
tanaka_psd,tanaka_dend = tanaka_data_s3()
# CF_PSD = get_CF_2023_4_data()
CF_PSD = pd.read_csv("./CF_2023_T4.csv")
# breakpoint()
CF_PSD["GluA2_Area_per_cluster"] = CF_PSD["GluA2_Area"]/CF_PSD["#_Cluster"]
CF_PSD["PSD_Area_per_cluster"] = CF_PSD["PSD_Area"]/CF_PSD["#_Cluster"]
CF_PSD.to_csv("./CF_2023_T4.csv", index=False)
CF_control_psd = CF_PSD.loc[CF_PSD['condition'] == "Control"]
CF_cLTP_psd = CF_PSD.loc[CF_PSD['condition'] == "cLTP"]
CF_time_points = CF_PSD["timepoint"].unique()
# breakpoint()

def plotStimuLocation(data_to_plot,op_folder,time_to_plot,locations,stim_start,stim_end,sou="s",ax_label=0):
    stim_or_unstim = "ST"
    if sou == "u":
        stim_or_unstim="UN"
    fig, ax = plt.subplots(figsize = (8,6),nrows=1,ncols=1)
    tps = time_to_plot/60
    ax.tick_params(axis='both', which='major', labelsize=18)
    # ax.set_ylim([75,145])
    for loc in locations:
        total_surf = data_to_plot[labels[1]][loc,:] #+data[labels[2]][loc, :]
        ax.plot(tps,100*(total_surf/total_surf[0]), color = COLORS_dict["shaft_s"], lw=3, label=r"$P_s$")
        ax.plot(tps, 100*(data_to_plot[labels[0]][loc, :]/data_to_plot[labels[0]][loc, 0]), color=COLORS_dict["shaft_i"], lw=3, label=r"$P_c$")
        ax.plot(tps, 100*(data_to_plot[labels[2]][loc, :]/data_to_plot[labels[2]][loc, 0]), color=COLORS_dict["spine_s"], lw=3, label=r"$P_{spine}$")
        if sou == "s":
            ax.fill_betweenx(range(85,125,1), stim_start/60,stim_end/60, alpha=0.2,color="r")

        ax.plot(tps,100*np.ones(tps.shape),'k--')
        x_loc = math.ceil(loc*dx)
        ax.set_xticks(np.arange(0,70,10))
        # ax.set_xticks(x_tics)
        if ax_label==1:
            ax.set_xlabel("Time (Minutes)")
            ax.set_ylabel("Normalized protein conc.")
        ax.spines[["right", "top"]].set_visible(False)
        ax.legend(loc="upper right", frameon=False, fontsize=f_size)
        # xlabels = [item.get_text() for item in ax.get_xticklabels()]
        # xlabels = [int(float(i)/60) for i in xlabels]
        # ax.set_xticklabels(xlabels)
    # plt.legend()
    # ax.errorbar(data_4g[:,0]+5,data_4g[:,1]/100,data_4g[:,2]/100,color = "orange")
        plt.tight_layout()
        SaveFigures("{}/{}_Time_Evolution_{}_at_location_{}".format(op_folder,stim_or_unstim,date_time,loc))
        plt.show()

plottinglabs = [r"$P_c$",r"$P_s$",r"$P_{spine}$"]
def plotIndividual(data1,data2,op_folder,tps1,tps2,lc1,lab_index,d1_suf,d2_suf,ax_label=0):
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    ax.plot(tps1,data1[labels[lab_index]][lc1, :]/data1[labels[lab_index]][lc1, 0],label=plottinglabs[lab_index]+d1_suf)
    ax.plot(tps2, data2[labels[lab_index]][lc1, :]/data2[labels[lab_index]][lc1, 0], label=plottinglabs[lab_index] + d2_suf)
    total = data1[labels[lab_index]][lc1, :] + data2[labels[lab_index]][lc1, :]
    ax.plot(tps2,total/total[0] ,
            label=plottinglabs[lab_index] + d1_suf+d2_suf)

    if ax_label == 1:
        ax.set_xlabel("Time (Minutes)")
        ax.set_ylabel("Normalized protein conc.")
    ax.spines[["right", "top"]].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=f_size)
    plt.tight_layout
    SaveFigures("{}/{}_Time_Evolution_{}_at_location_{}".format(op_folder, labels[lab_index], date_time, lc1))
    plt.show()

# plotIndividual(aligned_d1,aligned_d2,tp1,tp2,250,2," R1"," R2")
def comparePatterson(loc,op_folder):
    fig,ax = plt.subplots(figsize = (8,6),nrows=1,ncols=1)
    tps = (timepoints-30)/60
    total_surf = (data[labels[1]][loc,:] + data[labels[2]][loc, :])
    total_surf = total_surf[0]
    # breakpoint()
    perce_increas = (total_surf - total_surf[0])/total_surf[0]
    plt.errorbar(patterson_1g[:,0],patterson_1g[:,1],patterson_1g[:,2],color="k",marker='d',linestyle='--',label="data")
    plt.plot(tps,  total_surf, color="blue", label="model")
    # plt.ylabel(r"$\Delta$ Total Surf count ($\%$)")
    # plt.xlabel(r"Simulation Time (in mins)")
    ax.hlines(y=160, xmin=0, xmax=1, linewidth=2, color='r')
    ax.text(x=0, y=162, s="Stim", fontsize=f_size)
    plt.legend(fontsize=f_size)
    plt.tight_layout()
    SaveFigures("{}/Patterson_2010_fit_{}".format(op_folder, date_time))
    plt.show()


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = a.shape[0]
    m, se = np.mean(a,axis=1),sst.sem(a,axis=1)
    h = se * sst.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

def compareTanaka(data,fname):
    """
    Compares the plasticity induced changes in GLuA1 in synaptics and extrasynaptic compartments
    in model (namely, P_spine and P_s) and takes an average across the complete model dendrite
    agains the PLSM and non-PSLM changes reported in Tanaka et al., 2012 cell reports paper
    """
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    tps = (timepoints - 30) / 60
    # model_dend = (data[labels[1]][loc, :] / data[labels[1]][loc, 0])
    # model_psd = (data[labels[2]][loc, :] / data[labels[2]][loc, 0])
    model_dend, model_dend_cim,model_dend_cip = mean_confidence_interval(np.divide(data[labels[1]].T, data[labels[1]][:, 0]))
    model_psd,model_psd_cim,model_psd_cip = mean_confidence_interval(np.divide(data[labels[2]].T, data[labels[2]][:, 0]))

    ax.errorbar(tanaka_dend[:, 0], tanaka_dend[:, 1], tanaka_dend[:, 2], color="#6382f6fa", marker='^', linestyle='',
                 label="GluA1 in Extra syn.",markersize=20,linewidth=3,capsize=6,zorder=10)
    ax.errorbar(tanaka_psd[:, 0], tanaka_psd[:, 1], tanaka_psd[:, 2], color="#0135f8ff", marker='o', linestyle='',
                label="GluA1 in PSD", markersize=20,linewidth=3, capsize=6,zorder=10,)

    ax.plot(tps, 100 * model_psd, color="#0135f8ff", label=r"$P_{PSD}$: GluA1",linewidth=4)
    # ax.fill_between(tps, 100*model_dend_cim, 100*model_dend_cip,color=COLORS_dict["spine_s"], alpha=.3)
    ax.plot(tps, 100 * model_dend, color="#6382f6fa", label=r"$P_s$: GluA1",linewidth=4)
    # ax.fill_between(tps, 100 * model_psd_cim, 100 * model_psd_cip, color=COLORS_dict["shaft_s"], alpha=.3)
    # breakpoint()
    ax.hlines(y=100,xmin=-10,xmax=tps[-1],linestyles="dotted",color='k')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.axvspan(0, 10, alpha=0.5, color='gray')
    # ax.hlines(y=91, xmin=0, xmax=10, linewidth=2, color='r')
    # ax.text(x=0, y=93.0,s="Gly stim.",fontsize=f_size)#, transform=ax.transAxes)
    ax.set_xlabel("Time (min)",fontsize=f_size)
    ax.set_ylabel(r"$\%\Delta$ in GluA1 fluorescence  ",fontsize=f_size)
    ax.set_xticks(np.arange(-10,61,10))
    ax.set_ylim([90,150])
    ax.set_yticks(np.arange(90,151,10))
    ax.set_xlim([-12, 65])
    ax.tick_params(axis='both', which='major', labelsize=f_size)
    plt.legend(frameon=False,fontsize=f_size,loc="upper right",labelcolor='linecolor')
    plt.tight_layout()
    SaveFigures(fname)
    plt.show()

def compareClavet_Fournier(data,fname,init_time=30):
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    tps = (timepoints - init_time) / 60
    model_dend = np.divide(data[labels[1]].T, data[labels[1]][:, 0]).mean(axis=1)
    model_psd = np.divide(data[labels[2]].T, data[labels[2]][:, 0]).mean(axis=1)
    # breakpoint()
    ax.plot(tps, 100*model_psd, color="#017a56ff", label=r"$P_{PSD}$: GluA2")
    ax.plot(tps, 100*model_dend, color="#63bf1eff", label=r"$P_{s}: GluA2$")
    COI = "GluA2_Area"
    cLTP_means = CF_cLTP_psd.groupby('timepoint')[[COI]].mean().to_numpy()
    cLTP_stds = CF_cLTP_psd.groupby('timepoint')[[COI]].std().to_numpy()
    # cLTP_counts = CF_cLTP_psd.groupby('timepoint')[[COI]].count().to_numpy()

    cont_means = CF_control_psd.groupby('timepoint')[[COI]].mean().to_numpy()
    # breakpoint()
    cLTP_means /= cont_means
    cLTP_means /= cLTP_means[0]
    cont_means /= cont_means[0]
    # cont_stds = CF_control_psd.groupby('timepoint')[[COI]].std().to_numpy()
    # cont_counts = CF_control_psd.groupby('timepoint')[[COI]].count().to_numpy()
    # breakpoint()
    ax.errorbar(CF_time_points,100*cLTP_means[:,0],100*cLTP_stds[:,0], color="#017a56ff", marker='o', linestyle='',
                label="GluA2 in PSD",markersize=20,linewidth=3, capsize=6,zorder=10)
    # ax.errorbar(CF_time_points, cont_means[:, 0], cont_stds[:, 0], color="k", marker='^', linestyle='',
    #             label="{}: Cont".format(COI), markersize=12, capsize=3, zorder=10)
    ax.hlines(y=100,xmin=-10,xmax=tps[-1],linestyles="dotted",color='k')
    # ax.hlines(y=85, xmin=l_start / 60, xmax=l_end / 60, linewidth=5, color='k')
    # ax.text(x=l_start/60, y=90, s=r"$\uparrow$ "+mesg, fontsize=f_size)
    # ax.hlines(y=75, xmin=s_start/60, xmax=s_end/60, linewidth=5, color='r')
    # ax.text(x=0, y=72.0,s="Gly stim.",fontsize=f_size)
    # ax.axvspan(l_start/60, l_end/60, alpha=0.5, color='gray')
    ax.set_ylim([70,200])
    ax.set_xticks(np.arange(0,tps[-1]+1,60))
    ax.hlines(y=100, xmin=tps[0], xmax=tps[-1], linestyles="dotted", color='k')
    plt.legend(frameon=False,fontsize=f_size,loc="upper right")
    ax.set_xlabel("Time (min)", fontsize=f_size)
    ax.set_ylabel(r"$\%\Delta$ in GluA2 fluorescence  ", fontsize=f_size)
    ax.tick_params(axis='both', which='major', labelsize=f_size)
    plt.tight_layout()

    # print("op folder = ",op_folder)
    SaveFigures(fname)
    plt.show()
def Plot3DMatrixIndividual(dat,locn,off_set,step,stim_start,stim_end,title,lab,c_map,ax_label=0):
    min_y = locn - off_set
    max_y = locn + off_set
    x_g = [i for i in range(int(min_y / dx), int(max_y / dx), int(step / dx))]
    fps = 60
    num_frames = int(fps / dt)

    f_rames = [i for i in range(0, len(timepoints), num_frames)]
    pc_loc_data = dat[x_g, :]

    pc_init_data = pc_loc_data[:, 0]

    pc_time_loc_data = pc_loc_data[:, f_rames]

    pc_evol = (pc_time_loc_data.T / pc_init_data).T
    asp = 'auto'
    # f_size = 14
    xlab, ylab = "Time (min)", r"Distance from Soma ($\mu$m)"
    pos, orient, pad = 'right', 'vertical', 0.5
    fig, ax1 = plt.subplots(figsize=(8, 6), ncols=1)
    # ticks = np.arange(1, 2 * off_set + 10, 10)
    # ax1.set_yticks(ticks, np.arange(locn - off_set, locn + off_set + 10, 10))
    ax1.set_xticks(np.arange(0, 70,10))
    # ax1 = fig.add_subplot(131)
    # breakpoint()
    im1 = ax1.imshow(pc_evol, aspect=asp, cmap=c_map)
    if ax_label==1:
        ax1.set_ylabel(ylab, fontsize=f_size)
        ax1.set_xlabel(xlab, fontsize=f_size)
        ax1.set_title(title,fontsize=f_size)
    ax1.tick_params(labelsize=12)
    # ax1.set(yticks=y_positions, yticklabels=y_labels);
    # ax1.set(xticks=x_positions,xticklabels=f_rames)
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes(pos, size='5%', pad=pad)
    plt.colorbar(im1, cax=cax1, orientation=orient)
    ax1.hlines(y=20, xmin=stim_start / 60, xmax=stim_end / 60, linewidth=2, color='k')
    # ax1.text(x=5, y=10, s="Translation \n up-regulation", color="k", fontsize=f_size)
    # breakpoint()
    # locs, labs = ax1.yticks()
    # ax1.text(x=2, y=92.0, s="Gly stimulation", fontsize=16)  # , transform=ax.transAxes)
    # ax2 = fig.add_subplot(132)

    plt.tight_layout()
    SaveFigures("{}/Matrix_plot_{}_{}".format(op_folder,lab, date_time))
    plt.show()


def HeterosynapticAvgTemporal(data_to_plot,timepoints,stim_loc=[],unstim_loc=[]):
    tps = (timepoints-30)/60
    fig, ax = plt.subplots(figsize=(8, 6), ncols=1)
    mean_stim_ps = data_to_plot[labels[1]][stim_loc].mean(axis=0)  # +data[labels[2]][loc, :]
    mean_stim_pc = data_to_plot[labels[0]][stim_loc].mean(axis=0)
    mean_stim_pspine = data_to_plot[labels[2]][stim_loc].mean(axis=0)

    mean_unstim_ps = data_to_plot[labels[1]][unstim_loc].mean(axis=0)  # +data[labels[2]][loc, :]
    mean_unstim_pc = data_to_plot[labels[0]][unstim_loc].mean(axis=0)
    mean_unstim_pspine = data_to_plot[labels[2]][unstim_loc].mean(axis=0)
    # breakpoint()
    ax.plot(tps, 100 * (mean_stim_ps / mean_stim_ps[0]),
            color=COLORS_dict["shaft_s"], lw=3,linestyle='-', label=r"S:$P_s$")
    ax.plot(tps, 100 * (mean_stim_pc / mean_stim_pc[0]),
            color=COLORS_dict["shaft_i"], lw=3,linestyle='-', label=r"S:$P_c$")
    ax.plot(tps, 100 * (mean_stim_pspine/mean_stim_pspine[0]),
            color=COLORS_dict["spine_s"], lw=3,linestyle='-', label=r"S:$P_{spine}$")

    ax.plot(tps, 100 * (mean_unstim_ps / mean_unstim_ps[0]),
            color=COLORS_dict["shaft_s"], lw=3, linestyle='--', label=r"U:$P_s$")
    ax.plot(tps, 100 * (mean_unstim_pc / mean_unstim_pc[0]),
            color=COLORS_dict["shaft_i"], lw=3, linestyle='--', label=r"U:$P_c$")
    ax.plot(tps, 100 * (mean_unstim_pspine / mean_unstim_pspine[0]),
            color=COLORS_dict["spine_s"], lw=3, linestyle='--', label=r"U:$P_{spine}$")
    plt.plot(tps, 100*np.ones(tps.shape), color='black',linestyle='-.', lw=1, label=labels[0])
    ax.spines[["right", "top"]].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=f_size)
    plt.tight_layout()
    plt.show()

per = "100"
labels = ["PC", "PS", "PSPINE"]

if __name__ == '__main__':
    """
        parsing argumenst on mRNA and widths to be analysed. Each length is analysed seperately. mRNAs can be analysed in combinations
    """
    parser = argparse.ArgumentParser(description='mRNA analysis.py file ')
    parser.add_argument('-d', "--date_time", type=str, default="11_08_2024_16_18_55",
                        help='name of the folder to read data from')
    parser.add_argument('-s', "--subunit", type=str, default="GluA2",
                        help='name of the GluA subunit, GluA1 or GluA2')

    # reading the argumenst
    args = parser.parse_args()
    date_time = args.date_time
    subunit  = args.subunit
    # s_start = 0
    # s_end = 10*60
    # l_start = 0
    # l_end = 10*60

    root_folder = "./Time-dependent/GluAsims/GluA_simulations"
    file_names = ["{}_{}_{}_percent.npy".format(i, date_time, per) for i in
                  labels]  # ["PC_{}_{}_percent.npy".format(date_time,per),"PS_{}_{}_percent.npy".format(date_time, per),"PSPINE_{}_{}_percent.npy".format(date_time,per)]
    data = {}
    if subunit == "GluA2":
        input_folder = "{}/GluA2/CNIH2_dep/{}/".format(root_folder, date_time)
    else:
        input_folder = "{}/GluA1/{}/".format(root_folder, date_time)

    op_folder = os.path.join(input_folder + "figures_{}".format(date_time))
    os.makedirs(op_folder, exist_ok=True)

    for idx, fname in enumerate(file_names):
        print(labels[idx])
        data[labels[idx]] = np.load(input_folder + fname)
    timepoints = np.load(input_folder + "timepoints_{0}_{1}_percent.npy".format(date_time, per))
    dpi = 300
    interval = 500
    # locn = 50
    # beta_span = 3
    # alpha_span = 10
    # fps = 25
    f_size = 26
    P_s_init, P_c_init, P_spine_init, SP_model1 = RunModelWithFile(baseline_param_file)
    # L = 500.0
    dx = SP_model1.dx
    x_grid = SP_model1.x_grid  # np.arange(0,SP_model1L,dx)
    x_points = x_grid.shape[0]
    num_frames = int((0.25 * 60) / dt)
    f_rames = [i for i in range(0, len(timepoints), num_frames)]

    if subunit == "GluA1":
        compareTanaka(data,"{}/Tanaka_2012_fit_{}".format(op_folder, date_time))
    else:
        compareClavet_Fournier(data,"{}/Clavet_Fournier_2023_fit_{}".format(op_folder, date_time))
