#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:26:13 2022

@author: surbhitwagle
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

date_time = "01_16_2023_15_44_27"#"01_02_2023_13_13_25"#"12_25_2022_16_10_22"
per = "100"
labels = ["PC","PS","PSPINE"]
file_names = ["{}_{}_{}_percent.npy".format(i,date_time,per) for i in labels]#["PC_{}_{}_percent.npy".format(date_time,per),"PS_{}_{}_percent.npy".format(date_time, per),"PSPINE_{}_{}_percent.npy".format(date_time,per)]
input_folder = "{}/".format(date_time);
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
    
plt.xlabel('x'); 
plt.ylabel('concentration')

L = 500.0
dx = 0.24
x_grid = np.arange(0,L,dx)
x_points = x_grid.shape[0]

T = 100
dt = 1
t_grid =  np.arange(0,T+dt,dt)
t_points = t_grid.shape[0]

fig,ax = plt.subplots()
title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")
line, = ax.plot(x_grid,  data[labels[0]][0,:], color = 'blue', lw=1,label=labels[0])
line2, = ax.plot(x_grid, data[labels[1]][0,:], color = 'red', lw=1,label=labels[1])
line3, = ax.plot(x_grid, data[labels[2]][0,:], color = 'green', lw=1,label=labels[2])
plt.legend(loc="upper right")
f_rames = [i for i in range(0,t_points,1)]
def animate(i):
    line.set_data(x_grid,  data[labels[0]][i,:])
    line2.set_data(x_grid,  data[labels[1]][i,:])
    line3.set_data(x_grid,  data[labels[2]][i,:])
    title.set_text('time = %.3f s' % (i*dt))
    
    return  line, line2, line3, title
    # for k in range(t_points):
    #     for i,key in enumerate(labels):
    #         plots[i].set_data((x_grid, data[key][i,:]),label=key)
# breakpoint()
ani = FuncAnimation(fig, animate, interval=40, blit=True, repeat=True, frames=f_rames)    
ani.save("{}/{}.gif".format(input_folder,date_time), dpi=300, writer=PillowWriter(fps=25))
# plt.close()

fig,ax = plt.subplots()
line, = ax.plot(x_grid,  data[labels[0]][0,:]/data[labels[0]][0,:], color = 'blue', lw=1,label=labels[0])
line2, = ax.plot(x_grid, data[labels[1]][0,:]/data[labels[1]][0,:], color = 'red', lw=1,label=labels[1])
line3, = ax.plot(x_grid, data[labels[2]][0,:]/data[labels[2]][0,:], color = 'green', lw=1,label=labels[2])
plt.legend(loc="upper left")
plt.plot(x_grid,  np.ones(x_grid.shape), color = 'black', lw=1,label=labels[0])
# plt.ylim([0,1.5])
title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")
def animateRatio(i):
    line.set_data(x_grid,  data[labels[0]][i,:]/data[labels[0]][0,:])
    line2.set_data(x_grid,  data[labels[1]][i,:]/data[labels[1]][0,:])
    line3.set_data(x_grid,  data[labels[2]][i,:]/data[labels[2]][0,:])
    title.set_text('time = %.3f s'  % (i*dt))
    
    return  line, line2,line3,title

aniratio = FuncAnimation(fig, animateRatio, interval=40, blit=True, repeat=True, frames=f_rames)    
aniratio.save("{}/ratio_{}.gif".format(input_folder,date_time), dpi=300, writer=PillowWriter(fps=25))
# plt.close()
