#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 14:22:18 2021

@author: surbhitwagle
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

L = 221;    #lenght of the dendrite
class PlottingWidgetAMPA():
    def PlotSingleSimTwoProtein(self,x,ps,pc,lab_ps,lab_pc,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        ax.plot(x,ps,label=lab_ps)
        ax.plot(x,pc,label=lab_pc)
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        plt.title(title_string,fontsize=fsize)
        plt.legend()
        # plt.show()
        folder = "."
        if save_it == 1:
            plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=150)
            plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=150)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    def PlotSingleSimSingleProtein(self,x,y,lab,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        ax.plot(x,y,label=lab)
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        plt.title(title_string,fontsize=fsize)
        plt.legend()
        # plt.show()
        folder = "."
        if save_it == 1:
            plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=150)
            plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=150)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    # def PlotMultipleSim(self,x,y,lab,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
    #     # plotting multiple simulations
    #     # useful when trying single parameter manipulations
        
    #     if x.shape != y.shape:
    #         print("data shapes are not equal x is of size ", x.shape , "and y is of size ",y.shape )
    #     else:
    #         fig, ax = plt.subplots(figsize=(width, height))
    #         plt.rc('font', **{'family':'serif','serif':['Palatino']})
    #         plt.rc('text', usetex=True)
    #         y_shape = y.shape
    #         print(y_shape[0])
    #         for i in range(0,y_shape[0]):
    #             ax.plot(x[i],y[i],label=lab[i])
    #         ax.set_xlabel(xlab,fontsize=fsize)
    #         ax.set_ylabel(ylab,fontsize=fsize)
    #         plt.title(title_string,fontsize=fsize)
    #         plt.legend()
    #         # plt.show()
    #         folder = "."
    #         if save_it == 1:
    #             plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=150)
    #             plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=150)
    #             print("saved figures to: {%s/%s}" %(folder, file_name))
    #         else:
    #             print("Plots not saved")
    #         plt.show()