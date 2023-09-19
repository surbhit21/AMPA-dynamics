#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:26:13 2022

@author: surbhitwagle
"""
import math
import os.path

from AMPA_model import RunModelWithFile
from GetPublishedData import *

import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from Utility import *
from matplotlib.animation import FuncAnimation, PillowWriter

date_time = "09_18_2023_17_36_34"
per = "100"
labels = ["PC","PS","PSPINE"]
file_names = ["{}_{}_{}_percent.npy".format(i,date_time,per) for i in labels]#["PC_{}_{}_percent.npy".format(date_time,per),"PS_{}_{}_percent.npy".format(date_time, per),"PSPINE_{}_{}_percent.npy".format(date_time,per)]
input_folder = "./Time-dependent/{}/".format(date_time);
data = {}
op_folder = os.path.join(input_folder+"figures_{}".format(date_time))
os.makedirs(op_folder, exist_ok=True)
def SaveFigures(filename,ext_list = [".png",".svg",".pdf"],dpi=300):
    """

        function to save figures
        required arguments:
            filename
    """
    for ext in ext_list:
        plt.savefig(filename+ext,dpi=dpi)
            
            
for idx,fname in  enumerate(file_names):
    data[labels[idx]] = np.load(input_folder+fname)
timepoints = np.load(input_folder+"timepoints_{0}_{1}_percent.npy".format(date_time,per))
# plt.xlabel('x');
# plt.ylabel('concentration')
# breakpoint()
dpi = 300
interval = 500
locn = 50
beta_span = 3
alpha_span = 10
fps = 25
P_s_init,P_c_init,P_spine_init,SP_model1 = RunModelWithFile("./ModelParamsTemporal.json")
L = 500.0
dx = SP_model1.dx
x_grid = np.arange(0,L,dx)
x_points = x_grid.shape[0]
num_frames= int(2/0.02)
f_rames = [i for i in range(0,len(timepoints),num_frames)]
# T = 20
# dt = 0.002*10
# t_grid =  np.arange(0,T+dt,dt)
# t_points = t_grid.shape[0]
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
# txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# plt.legend(loc="upper right")
# plt.ylim([0,60])
# def animate(i):
#     line.set_data(x_grid,  data[labels[0]][:,i])
#     line2.set_data(x_grid,  data[labels[1]][:,i])
#     line3.set_data(x_grid,  data[labels[2]][:,i])
#     title.set_text('time = %.3f mins'  % (timepoints[i]/60))
#     txt1.set_text(r"stim location = {} $\mu m$".format(locn))
#     return  line, line2, line3, title,txt1
#
# ani = FuncAnimation(fig, animate, interval=interval, blit=True, repeat=True, frames=f_rames)
# ani.save("{}/{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()

"""
Plotting a time gif of ratio of change from initial state
"""

# fig,ax = plt.subplots()
# line, = ax.plot(x_grid,  data[labels[0]][:,0]/P_c_init, color = COLORS_dict["shaft_i"], lw=1,label=labels[0])
# line2, = ax.plot(x_grid, data[labels[1]][:,0]/P_s_init, color = COLORS_dict["shaft_s"], lw=1,label=labels[1])
# line3, = ax.plot(x_grid, data[labels[2]][:,0]/P_spine_init, color = COLORS_dict["spine_s"], lw=1,label=labels[2])
# txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# plt.legend(loc="upper left")
# plt.plot(x_grid,  np.ones(x_grid.shape), color = 'black', lw=1,label=labels[0])
# plt.ylim([0,2])
# title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# def animateRatio(i):
#     line.set_data(x_grid,  data[labels[0]][:,i]/P_c_init)
#     line2.set_data(x_grid,  data[labels[1]][:,i]/P_s_init)
#     line3.set_data(x_grid,  data[labels[2]][:,i]/P_spine_init)
#     title.set_text('time = %.3f mins'  % (timepoints[i]/60))
#     txt1.set_text(r"stim location = {} $\mu m$".format(locn))
#     return line, line2, line3, title, txt1
#
# aniratio = FuncAnimation(fig, animateRatio, interval=interval, blit=True, repeat=True, frames=f_rames)
# aniratio.save("{}/ratio_{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()


"""
Plotting a time gif of zoomed at location of stimlation loc with a surrounding area of span 
"""
# fig,ax = plt.subplots()
# # loc =  50
# span = 15
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
def plotStimuLocation(locations,sou="s"):
    stim_or_unstim = "ST"
    if sou == "u":
        stim_or_unstim="UN"
    fig, ax = plt.subplots(figsize = (8,6),nrows=1,ncols=1)
    tps = timepoints/60
    for loc in locations:
        total_surf = data[labels[1]][loc,:] #+data[labels[2]][loc, :]
        ax.plot(tps,100*(total_surf/total_surf[0]), color = COLORS_dict["shaft_s"], lw=1, label=r"$P_s$")
        ax.plot(tps, 100*(data[labels[0]][loc, :]/data[labels[0]][loc, 0]), color=COLORS_dict["shaft_i"], lw=1, label=r"$P_c$")
        ax.plot(tps, 100*(data[labels[2]][loc, :]/data[labels[2]][loc, 0]), color=COLORS_dict["spine_s"], lw=1, label=r"$P_{spine}$")
        if sou == "s":
            ax.fill_betweenx(range(85,125,1), 30/60,90/60, alpha=0.2,color="r")

        ax.plot(tps,100*np.ones(tps.shape),'k--')
        x_loc = math.ceil(loc*dx)
        ax.spines[["right","top"]].set_visible(False)
        ax.text(x=0.3,y=0.8,s="location = {}".format(x_loc),transform=ax.transAxes)
    # for axs in ax:
    #     ax.set_xlabel("T (Minutes)")
    #     ax.set_ylabel("Normalized Protein count")
        ax.legend(loc="upper right",frameon=False)
        # xlabels = [item.get_text() for item in ax.get_xticklabels()]
        # xlabels = [int(float(i)/60) for i in xlabels]
        # ax.set_xticklabels(xlabels)
    # plt.legend()
    # ax.errorbar(data_4g[:,0]+5,data_4g[:,1]/100,data_4g[:,2]/100,color = "orange")
    plt.tight_layout()
    SaveFigures("{}/{}_Time_Evolution_{}_at_location_{}".format(op_folder,stim_or_unstim,date_time,x_loc))
    plt.show()

def comparePatterson(loc):
    fig,ax = plt.subplots(figsize = (8,6),nrows=1,ncols=1)
    tps = (timepoints-30)/60
    total_surf = (data[labels[1]][loc,:] + data[labels[2]][loc, :])
    total_surf = total_surf[0]
    breakpoint()
    perce_increas = (total_surf - total_surf[0])/total_surf[0]
    plt.errorbar(patterson_1g[:,0],patterson_1g[:,1],patterson_1g[:,2],color="k",marker='d',linestyle='--',label="data")
    plt.plot(tps,  total_surf, color="blue", label="model")
    plt.ylabel(r"$\Delta$ Total Surf count ($\%$)")
    plt.xlabel(r"Simulation Time (in mins)")
    ax.hlines(y=160, xmin=0, xmax=1, linewidth=2, color='r')
    ax.text(x=0, y=162, s="Stim", fontsize=16)
    plt.legend(fontsize=16)
    plt.tight_layout()
    SaveFigures("{}/Patterson_2010_fit_{}".format(op_folder, date_time))
    plt.show()

def compareTanaka(loc):
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    tps = (timepoints - 30) / 60
    model_dend = (data[labels[1]][loc, :] / data[labels[1]][loc, 0])
    model_psd = (data[labels[2]][loc, :] / data[labels[2]][loc, 0])
    ax.plot(tps, 100 * model_psd[0], color=COLORS_dict["spine_s"], label=r"model:$P_{spine}$")
    ax.plot(tps, 100 * model_dend[0], color=COLORS_dict["shaft_s"], label=r"model:$P_s$", alpha=0.5)
    ax.errorbar(tanaka_psd[:, 0], tanaka_psd[:, 1], tanaka_psd[:, 2], color="k", marker='d', linestyle='--',
                 label="data:PSLM")
    ax.errorbar(tanaka_dend[:, 0], tanaka_dend[:, 1], tanaka_dend[:, 2], color="k", marker='d', linestyle='--',
                 label="data:N-PSLM")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.hlines(y=90, xmin=0, xmax=10, linewidth=2, color='r')
    ax.text(x=2, y=92.0,s="Gly stimulation",fontsize=16)#, transform=ax.transAxes)
    # ax.set_xlabel("Time in minutes")
    # ax.set_ylabel("Normalized protein count")

    plt.legend(frameon=False)
    plt.tight_layout()
    SaveFigures("{}/Tanaka_2012_fit_{}".format(op_folder, date_time))
    plt.show()

def Plot3DMatrix(locn,off_set,step,stim_start,stim_end):
    min_y = locn-off_set
    max_y = locn+off_set
    step_y = step
    x_g= [i for i in range(int(min_y/dx),int(max_y/dx),int(dx/dx))]
    fps = 5
    num_frames = int(fps / 0.02)
    stim_start /= fps
    stim_start/= num_frames
    stim_end /= fps
    stim_end /= num_frames

    f_rames = [i for i in range(0, len(timepoints), num_frames)]

    x = np.arange(0,timepoints[-1]+1,num_frames)
    nx = x.shape[0]
    no_labelsx = int((timepoints[-1] - 0) / num_frames)  # how many labels to see on axis x
    s_x = int(nx / (no_labelsx - 1))
    x_positions = fps*np.arange(0,x.shape[0],s_x)
    # breakpoint()
    y = np.arange(min_y, max_y+1, step_y)  # the grid to which your data corresponds
    ny = y.shape[0]
    no_labels = int((max_y-min_y)/step_y)  # how many labels to see on axis x
    s_y = int(ny / (no_labels - 1))  # step between consecutive labels
    y_positions = step_y*np.arange(0, ny, s_y)  # pixel count at label position
    y_labels = y # labels you want to see
    # plt.xticks(x_positions, x_labels)
    pc_loc_data = data[labels[0]][x_g,:]
    ps_loc_data = data[labels[1]][x_g, :]
    psd_loc_data = data[labels[2]][x_g, :]


    pc_init_data = pc_loc_data[:,0]
    ps_init_data = ps_loc_data[:,0]
    psd_init_data = psd_loc_data[:,0]


    pc_time_loc_data = pc_loc_data[:,f_rames]
    ps_time_loc_data = ps_loc_data[:,f_rames]
    psd_time_loc_data = psd_loc_data[:,f_rames]

    pc_evol = (pc_time_loc_data.T / pc_init_data).T
    ps_evol =  (ps_time_loc_data.T / ps_init_data).T
    psd_evol =  (psd_time_loc_data.T / psd_init_data).T
    asp = 1
    f_size = 14
    xlab,ylab = "Time in sec",r"Distance from Soma $in \mu m$"
    pos,orient,pad = 'bottom','horizontal',0.5
    fig,(ax1,ax2,ax3) = plt.subplots(figsize=(18, 6),ncols=3)
    # ax1 = fig.add_subplot(131)
    # breakpoint()
    im1 = ax1.imshow(pc_evol,aspect=asp,cmap='autumn')
    ax1.set_ylabel(ylab,fontsize=f_size)
    ax1.set_xlabel(xlab,fontsize=f_size)
    ax1.set_title(r"$P_c$")
    # ax1.set(yticks=y_positions, yticklabels=y_labels);
    # ax1.set(xticks=x_positions,xticklabels=f_rames)
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes(pos, size='5%', pad=pad)
    plt.colorbar(im1, cax=cax1, orientation=orient)
    ax1.hlines(y=1.1, xmin=stim_start, xmax=stim_end, linewidth=2, color='k',transform=ax1.transAxes)
    # breakpoint()
    # locs, labs = ax1.yticks()
    # ax1.text(x=2, y=92.0, s="Gly stimulation", fontsize=16)  # , transform=ax.transAxes)
    # ax2 = fig.add_subplot(132)
    im2 = ax2.imshow(ps_evol, interpolation='None',aspect=asp,cmap="Blues")
    ax2.set_ylabel(ylab,fontsize=f_size)
    ax2.set_xlabel(xlab,fontsize=f_size)
    ax2.set_title(r"$P_s$")
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes(pos, size='5%', pad=pad)
    plt.colorbar(im2, cax=cax2, orientation=orient);
    # ax2.set(yticks=y_positions, yticklabels=y_labels);
    # ax2.set(xticks=x_positions, xticklabels=f_rames)
    ax2.hlines(y=1.1,  xmin=stim_start, xmax=stim_end, linewidth=2, color='k',transform=ax2.transAxes)

    # ax2.text(x=2, y=92.0, s="Gly stimulation", fontsize=16)  # , transform=ax.transAxes)
    # ax3 = fig.add_subplot(133)

    im3 = ax3.imshow(psd_evol, interpolation='None',aspect=asp,cmap="RdPu")
    ax3.set_ylabel(ylab,fontsize=f_size)
    ax3.set_xlabel(xlab,fontsize=f_size)
    ax3.set_title(r"$P_{spine}$")
    # ax3.set(yticks=y_positions, yticklabels=y_labels);
    # ax3.set(xticks=x_positions, xticklabels=f_rames)
    divider3 = make_axes_locatable(ax3)
    cax3 = divider3.append_axes('bottom', size='5%', pad=.5)
    ax3.hlines(y=1.1, xmin=stim_start, xmax=stim_end, linewidth=2, color='k',transform=ax3.transAxes)
    # ax3.text(x=2, y=92.0, s="Gly stimulation", fontsize=16)  # , transform=ax.transAxes)
    plt.colorbar(im3, cax=cax3, orientation='horizontal');
    plt.tight_layout()
    SaveFigures("{}/Matrix_plot_{}".format(op_folder, date_time))
    plt.show()
    # breakpoint()
locs = [50]#i for i in range(250,266,1)]
for lo in locs:
    plotStimuLocation([int(lo/dx)],"s")
un_stim_loc = [45]
for lo in un_stim_loc:
    plotStimuLocation([int(lo/dx)],"u")
# plotStimuLocation([int(251/dx)])
# plotStimuLocation([int(240/dx)])
# plotStimuLocation([int(262/dx)])


# # plotStimuLocation([int(locn/dx)])
# plotStimuLocation([int((locn-3)/dx)])
# plotStimuLocation([int(46/dx)])
# plotStimuLocation([int(58/dx)])
# comparePatterson([int(locn/dx)])
# breakpoint()
# compareTanaka([int(locn/dx)])
Plot3DMatrix(locn,50,dx,30,150)

# def animateflux()