#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:26:13 2022

@author: surbhitwagle
"""
from AMPA_model import RunModelWithFile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

date_time = "08_17_2023_15_20_41"#"01_02_2023_13_13_25"#"12_25_2022_16_10_22"
per = "100"
labels = ["PC","PS","PSPINE"]
file_names = ["{}_{}_{}_percent.npy".format(i,date_time,per) for i in labels]#["PC_{}_{}_percent.npy".format(date_time,per),"PS_{}_{}_percent.npy".format(date_time, per),"PSPINE_{}_{}_percent.npy".format(date_time,per)]
input_folder = "./Time-dependent/{}/".format(date_time);
data = {}

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
interval = 100
fps = 25
P_s_init,P_c_init,P_spine_init,SP_model1 = RunModelWithFile("./ModelParams.json")
L = 500.0
dx = SP_model1.dx
x_grid = np.arange(0,L,dx)
x_points = x_grid.shape[0]
num_frames= int(10/0.002)
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
# plt.legend(loc="upper right")
# plt.ylim([0,160])
# def animate(i):
#     line.set_data(x_grid,  data[labels[0]][:,i])
#     line2.set_data(x_grid,  data[labels[1]][:,i])
#     line3.set_data(x_grid,  data[labels[2]][:,i])
#     title.set_text('time = %.3f mins'  % (timepoints[i]/60))
#
#     return  line, line2, line3, title
#
# ani = FuncAnimation(fig, animate, interval=interval, blit=True, repeat=True, frames=f_rames)
# ani.save("{}/{}.gif".format(input_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()

"""
Plotting a time gif of ratio of change from initial state
"""
# fig,ax = plt.subplots()
# line, = ax.plot(x_grid,  data[labels[0]][:,0]/P_c_init, color = 'blue', lw=1,label=labels[0])
# line2, = ax.plot(x_grid, data[labels[1]][:,0]/P_s_init, color = 'red', lw=1,label=labels[1])
# line3, = ax.plot(x_grid, data[labels[2]][:,0]/P_spine_init, color = 'green', lw=1,label=labels[2])
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
#
#     return  line, line2,line3,title
#
# aniratio = FuncAnimation(fig, animateRatio, interval=interval, blit=True, repeat=True, frames=f_rames)
# aniratio.save("{}/ratio_{}.gif".format(input_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()


"""
Plotting a time gif of zoomed at location of stimlation loc with a surrounding area of span 
"""
# fig,ax = plt.subplots()
# loc = 100
# span = 10
# index_min = int((loc - span)/dx)
# index_max = int((loc + span)/dx)
# # range_new = [index_min:index_max]
# x_grid_new = x_grid[index_min:index_max]
# title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# line, = ax.plot(x_grid_new,  data[labels[0]][index_min:index_max,0], color = 'blue', lw=1,label=labels[0])
# line2, = ax.plot(x_grid_new, data[labels[1]][index_min:index_max,0], color = 'red', lw=1,label=labels[1])
# line3, = ax.plot(x_grid_new, data[labels[2]][index_min:index_max,0], color = 'green', lw=1,label=labels[2])
# plt.legend(loc="upper right")
# plt.ylim([0,160])
#
# def animateZoomed(i):
#     line.set_data(x_grid_new, data[labels[0]][index_min:index_max, i])
#     line2.set_data(x_grid_new, data[labels[1]][index_min:index_max, i])
#     line3.set_data(x_grid_new, data[labels[2]][index_min:index_max, i])
#     title.set_text('time = %.3f mins' % (timepoints[i]/60))
#     return line,line2,line3,title
#
# anizoomed = FuncAnimation(fig, animateZoomed, interval=interval, blit=True, repeat=True, frames=f_rames)
# anizoomed.save("{}/zoomed_{}.gif".format(input_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()

"""
Plotting the time evolution of stimulated locations over time
"""
# print(timepoints.shape)
# breakpoint()
def plotStimuLocation(locations):
    fig, ax = plt.subplots(figsize = (8,6),nrows=1,ncols=1)
    tps = timepoints/60
    for loc in locations:
        ax.plot(tps,data[labels[0]][loc,:]/data[labels[0]][loc,0], color = 'blue', lw=1,label=labels[0])
        ax.plot(tps, data[labels[1]][loc, :]/data[labels[1]][loc,0], color='red', lw=1, label=labels[1])
        ax.plot(tps, data[labels[2]][loc, :]/data[labels[2]][loc,0], color='green', lw=1, label=labels[2])
    # for axs in ax:
        ax.set_xlabel("T (Minutes)")
        ax.set_ylabel("Normalized Protein count")
        ax.legend()
        # xlabels = [item.get_text() for item in ax.get_xticklabels()]
        # xlabels = [int(float(i)/60) for i in xlabels]
        # ax.set_xticklabels(xlabels)
    # plt.legend()
    plt.tight_layout()
    SaveFigures("{}/Time_Evolution_{}".format(input_folder,date_time))
    plt.show()
plotStimuLocation([int(100/dx)])