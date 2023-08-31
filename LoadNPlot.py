#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:26:13 2022

@author: surbhitwagle
"""
import os.path

from AMPA_model import RunModelWithFile
from GetPublishedData import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

date_time = "08_31_2023_14_45_34"
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
plasticity_span = 3
fps = 25
P_s_init,P_c_init,P_spine_init,SP_model1 = RunModelWithFile("./ModelParams.json")
L = 500.0
dx = SP_model1.dx
x_grid = np.arange(0,L,dx)
x_points = x_grid.shape[0]
num_frames= int(5/0.02)
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
# line, = ax.plot(x_grid,  data[labels[0]][:,0], color = 'blue', lw=1,label=labels[0])
# line2, = ax.plot(x_grid, data[labels[1]][:,0], color = 'red', lw=1,label=labels[1])
# line3, = ax.plot(x_grid, data[labels[2]][:,0], color = 'green', lw=1,label=labels[2])
#
# txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# plt.legend(loc="upper right")
# plt.ylim([0,160])
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
# line, = ax.plot(x_grid,  data[labels[0]][:,0]/P_c_init, color = 'blue', lw=1,label=labels[0])
# line2, = ax.plot(x_grid, data[labels[1]][:,0]/P_s_init, color = 'red', lw=1,label=labels[1])
# line3, = ax.plot(x_grid, data[labels[2]][:,0]/P_spine_init, color = 'green', lw=1,label=labels[2])
# txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# plt.legend(loc="upper left")
# plt.plot(x_grid,  np.ones(x_grid.shape), color = 'black', lw=1,label=labels[0])
# plt.ylim([0,5])
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
fig,ax = plt.subplots()
# loc =  50
span = 10
index_min = int((locn - span)/dx)
index_max = int((locn + span)/dx)
plasticity_minx = int((locn-plasticity_span)/dx)
plasticity_maxx = int((locn+plasticity_span)/dx)
# range_new = [index_min:index_max]
x_grid_new = x_grid[index_min:index_max]
title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")
line, = ax.plot(x_grid_new,  data[labels[0]][index_min:index_max,0], color = 'blue', lw=1,label=labels[0])
line2, = ax.plot(x_grid_new, data[labels[1]][index_min:index_max,0], color = 'red', lw=1,label=labels[1])
line3, = ax.plot(x_grid_new, data[labels[2]][index_min:index_max,0], color = 'green', lw=1,label=labels[2])
ax.plot(x_grid_new,np.ones(x_grid_new.shape),"k--")
txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
min_line = ax.vlines(x= x_grid[plasticity_minx],ymin=0.5,ymax=1.4,linestyles="dashed")
max_line = ax.vlines(x= x_grid[plasticity_maxx],ymin=0.5,ymax=1.4,linestyles="dashed")
plt.legend(loc="upper right")
plt.ylim([0.5,1.4])

def animateZoomed(i):
    line.set_data(x_grid_new, data[labels[0]][index_min:index_max, i]/data[labels[0]][index_min:index_max, 0])
    line2.set_data(x_grid_new, data[labels[1]][index_min:index_max, i]/data[labels[1]][index_min:index_max, 0])
    line3.set_data(x_grid_new, data[labels[2]][index_min:index_max, i]/data[labels[2]][index_min:index_max, 0])
    title.set_text('time = %.3f mins' % ((timepoints[i]-30)/60))
    txt1.set_text(r"stim location = {} $\mu m$".format(locn))
    return line, line2, line3, title, txt1

anizoomed = FuncAnimation(fig, animateZoomed, interval=interval, blit=True, repeat=True, frames=f_rames)
anizoomed.save("{}/zoomed_{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
plt.close()

"""
Plotting the time evolution of stimulated locations over time
"""
# print(timepoints.shape)
# breakpoint()
patterson_1g = patterson_data_1g()
patterson_4c = patterson_data_4c()
tanaka_psd,tanaka_dend = tanaka_data_s3()
def plotStimuLocation(locations):
    fig, ax = plt.subplots(figsize = (8,6),nrows=1,ncols=1)
    tps = timepoints/60
    for loc in locations:
        total_surf = data[labels[1]][loc,:] #+data[labels[2]][loc, :]
        ax.plot(tps,100*(total_surf/total_surf[0]), color = 'blue', lw=1,label=labels[1])
        ax.plot(tps, 100*(data[labels[0]][loc, :]/data[labels[0]][loc, 0]), color='red', lw=1, label=labels[0])
        ax.plot(tps, 100*(data[labels[2]][loc, :]/data[labels[2]][loc, 0]), color='green', lw=1, label=labels[2])
        ax.plot(tps,100*np.ones(tps.shape),'k--')
        ax.text(x=0.3,y=0.8,s="location = {}".format(locn),transform=ax.transAxes)
    # for axs in ax:
        ax.set_xlabel("T (Minutes)")
        ax.set_ylabel("Normalized Protein count")
        ax.legend()
        # xlabels = [item.get_text() for item in ax.get_xticklabels()]
        # xlabels = [int(float(i)/60) for i in xlabels]
        # ax.set_xticklabels(xlabels)
    # plt.legend()
    # ax.errorbar(data_4g[:,0]+5,data_4g[:,1]/100,data_4g[:,2]/100,color = "orange")
    plt.tight_layout()
    SaveFigures("{}/Time_Evolution_{}".format(op_folder,date_time))
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
    ax.plot(tps, 100 * model_psd[0], color="blue", label=r"model:$P_{spine}$")
    ax.plot(tps, 100 * model_dend[0], color="blue", label=r"model:$P_s$", alpha=0.5)
    ax.errorbar(tanaka_psd[:, 0], tanaka_psd[:, 1], tanaka_psd[:, 2], color="k", marker='d', linestyle='--',
                 label="data:PSLM")
    ax.errorbar(tanaka_dend[:, 0], tanaka_dend[:, 1], tanaka_dend[:, 2], color="k", marker='d', linestyle='--',
                 label="data:N-PSLM",alpha=0.5)

    ax.hlines(y=90, xmin=0, xmax=10, linewidth=2, color='r')
    ax.text(x=2, y=92.0,s="Gly stimulation",fontsize=16)#, transform=ax.transAxes)
    ax.set_xlabel("Time in minutes")
    ax.set_ylabel("Normalized protein count")

    plt.legend()
    plt.tight_layout()
    SaveFigures("{}/Tanaka_2012_fit_{}".format(op_folder, date_time))
    plt.show()

plotStimuLocation([int(250/dx)])
# comparePatterson([int(locn/dx)])
# breakpoint()
compareTanaka([int(250/dx)])

# def animateflux()