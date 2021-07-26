#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 14:22:18 2021

@author: surbhitwagle
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
from pathlib import Path
L = 221;    #lenght of the dendrite
CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
color_list = [CB91_Blue, CB91_Pink, CB91_Green, CB91_Amber, CB91_Purple, CB91_Violet]

class PlottingWidgetAMPA():
    def CreateFolderRecursive(self,folder):
        Path(folder).mkdir(parents=True, exist_ok=True)
    def PlotSingleSimTwoProtein(self,x,ps,pc,lab_ps,lab_pc,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        ax.plot(x,ps,label=lab_ps,color=color_list[0])
        ax.plot(x,pc,label=lab_pc,color=color_list[1])
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
    
    def PlotMultiSimTwoProtein(self,x,ps,pc,lab_ps,lab_pc,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
        if x.shape != (ps[0].shape) and (ps[0].shape) != (pc[0].shape):
            print("data shapes are not equal x is of size ")
        else:
            markers = ['--r','-b',':g','_c','-.y'];
            fig, ax = plt.subplots(figsize=(width, height))
            # plt.rc('font', **{'family':'serif','serif':['Palatino']})
            # plt.rc('text', usetex=False)
            y_shape = ps.shape;
            for i in range(0,y_shape[0]):
                ax.plot(x,ps[i],markers[i],label=lab_ps[i],alpha=0.9,color=color_list[i])
                ax.plot(x,pc[i],markers[i],alpha=0.2,color=color_list[i])
            ax.set_xlabel(xlab,fontsize=fsize)
            ax.set_ylabel(ylab,fontsize=fsize)
            plt.title(title_string,fontsize=fsize)
            plt.legend()
            # plt.show()
            folder = "."
            print(file_name)
            if save_it == 1:
                plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=150)
                plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=150)
                print("saved figures to: {%s/%s}" %(folder, file_name))
            else:
                print("Plots not saved")
            plt.show()
    
    def PlotSingleSimSingleProtein(self,x,y,lab,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
        print(width,height)
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        ax.plot(x,y,label=lab,color=color_list[0])
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
    
    def PlotMultipleSim(self,x,y,lab=None,xlab=None,ylab=None,title_string=None,file_name=None,width=8,height=8,fsize=16,save_it = 1):
        # plotting multiple simulations
        # useful when trying single parameter manipulations
        
        if x.shape != (y[0].shape):
            print("data shapes are not equal x is of size ", x.shape , "and y is of size ",y.shape )
        else:
            fig, ax = plt.subplots(figsize=(width, height))
            plt.rc('font', **{'family':'serif','serif':['Palatino']})
            plt.rc('text', usetex=True)
            y_shape = y.shape
            print(y_shape[0])
            for i in range(0,y_shape[0]):
                ax.plot(x,y[i],label=lab[i],color=color_list[i])
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
        
    def plotNormSingle(self,x,y,lab,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1,norm_type = "Sum"):
        norm_y = y;
        new_ylab = ylab;
        if norm_type == "Sum":
            norm_y = y/np.sum(y);
            new_ylab = "Sum normalized " + ylab;
        elif norm_type == "Max":
            norm_y = y/np.max(y);
            new_ylab = "Max normalized " + ylab;
        elif norm_type == "Min":
            norm_y = y/np.min(y);
            new_ylab = "Min normalized " + ylab;
        elif norm_type == "Mean":
            norm_y = y/np.mean(y);
            new_ylab = "Mean normalized " + ylab;
        elif norm_type == "First":
            norm_y = y/y[0];
            new_ylab = "First value normalized " + ylab;
        elif norm_type == "Last":
            norm_y = y/y[-1];
            new_ylab = "Last value normalized " + ylab;
        new_file_name = file_name+"_"+norm_type;
        # print(new_file_name)
        # print(width,height)
        self.PlotSingleSimSingleProtein(x,norm_y,lab,xlab,new_ylab,title_string,new_file_name,width,height,fsize,save_it)
    
    def plotNormTwo(self,x,y1,y2,lab1,lab2,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1,norm_type = "Sum"):
        norm_y1 = y1;
        norm_y2 = y2;
        new_ylab = ylab;
        if norm_type == "Sum":
            norm_y1 = y1/np.sum(y1);
            norm_y2 = y2/np.sum(y1);
            new_ylab = "Sum normalized " + ylab;
        elif norm_type == "Max":
            norm_y1 = y1/np.max(y1);
            norm_y2 = y2/np.max(y1);
            new_ylab = "Max normalized " + ylab;
        elif norm_type == "Min":
            norm_y1 = y1/np.min(y1);
            norm_y2 = y2/np.min(y1);
            new_ylab = "Min normalized " + ylab;
        elif norm_type == "Mean":
            norm_y1 = y1/np.mean(y1);
            norm_y2 = y2/np.mean(y1);
            new_ylab = "Mean normalized " + ylab;
        elif norm_type == "First":
            norm_y1 = y1/y1[0];
            norm_y2 = y2/y1[0];
            new_ylab = "First value normalized " + ylab;
        elif norm_type == "Last":
            norm_y1 = y1/y1[-1];
            norm_y2 = y2/y1[-1];
            new_ylab = "Last value normalized " + ylab;
        new_file_name = file_name+"_"+norm_type;
        # print(new_file_name)
        # print(width,height)
        self.PlotSingleSimTwoProtein(x,norm_y1,norm_y2,lab1,lab2,xlab,new_ylab,title_string,new_file_name,width,height,fsize,save_it)
        
    def plotNormMulti(self,x,y,lab,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1,norm_type = "Sum"):
        norm_factor = 1;
        new_ylab = ylab;
        if norm_type == "Sum":
             norm_factor = y.sum(axis=1)
             new_ylab = "Sum normalized " + ylab;
        elif norm_type == "Max":
             norm_factor= y.max(axis=1)
             new_ylab = "Max normalized " + ylab;
        elif norm_type == "Min":
            norm_factor = y.min(axis=1);
            new_ylab = "Min normalized " + ylab;
        elif norm_type == "Mean":
            new_ylab = "Mean normalized " + ylab;
        elif norm_type == "First":
            norm_factor = y[:,0];
            new_ylab = "First value normalized " + ylab;
        elif norm_type == "Last":
            norm_factor = y[:,-1];
            new_ylab = "Last value normalized " + ylab;
        norm_y = y/norm_factor[:, np.newaxis]
        new_file_name = file_name+"_"+norm_type;
        # print(new_file_name)
        # print(width,height)
        self.PlotMultipleSim(x,norm_y,lab,xlab,new_ylab,title_string,new_file_name,width,height,fsize,save_it)
    
    def plotNormMultiTwo(self,x,y1,y2,lab1,lab2,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1,norm_type = "Sum"):
        norm_factor = 1;
        new_ylab = ylab;
        if norm_type == "Sum":
             norm_factor = y1.sum(axis=1)
             new_ylab = "Sum normalized " + ylab;
        elif norm_type == "Max":
             norm_factor= y1.max(axis=1)
             new_ylab = "Max normalized " + ylab;
        elif norm_type == "Min":
            norm_factor = y1.min(axis=1);
            new_ylab = "Min normalized " + ylab;
        elif norm_type == "Mean":
            new_ylab = "Mean normalized " + ylab;
        elif norm_type == "First":
            norm_factor = y1[:,0];
            new_ylab = "First value normalized " + ylab;
        elif norm_type == "Last":
            norm_factor = y1[:,-1];
            new_ylab = "Last value normalized " + ylab;
        norm_y1 = y1/norm_factor[:, np.newaxis]
        norm_y2 = y2/norm_factor[:, np.newaxis]
        new_file_name = file_name+"_"+norm_type;
        # print(new_file_name)
        # print(width,height)
        self.PlotMultiSimTwoProtein(x,norm_y1,norm_y2,lab1,lab2,xlab,new_ylab,title_string,new_file_name,width,height,fsize,save_it)