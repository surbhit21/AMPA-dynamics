#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 17:10:24 2022

@author: surbhitwagle
"""

# from AMPA_model import *
import csv
from functools import reduce
import json
from lmfit import conf_interval, minimize,Minimizer, Parameters, Parameter, report_fit, printfuncs
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import os
# import PlotBinned 
# from PlottingWidgetAMPA import *
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
from pathlib import Path

from scipy.optimize import curve_fit
from scipy.stats import ks_2samp, kruskal

# import scikit_posthocs as sps
import seaborn as sns
import scikit_posthocs as sp

COLORS = ["#005f73","#9b2226","#CA6702","#337357"]

class SpineDendriteRatio():
    def __init__(self,DataFolder):
        self.df = DataFolder
    def Readfiles(self,exclude_cells = []):
        in_folder = self.df+"merged/"
        files = os.listdir(in_folder) 
        # print(files)
        dend_data = {}
        spine_data = {}
        for f in files:
            if f.startswith('cell_'):
                f_dend = open(in_folder+f+"/dend_stat.json")
                dend_data[f] = json.load(f_dend)
                f_dend.close()
                f_spine = open(in_folder+f+"/Synapse_l.json")
                spine_data[f] = json.load(f_spine)
                f_spine.close()
        return dend_data,spine_data
    def GetRatios(self,measure,dend_data,spine_data):
        dend_rid = {}
        spine_rid = {}
        ratio_rid = []
        measure = 'RawIntDen'
        for cell in dend_data.keys():
            dend_rid[cell] = np.array(dend_data[cell]["Dendrite0"][measure])
            dend_rid[cell] = dend_rid[cell]/dend_rid[cell][0]
            spine_rid[cell] = np.zeros(dend_rid[cell].shape)
            for idx in range(len(spine_data[cell])):
                spine_rid[cell] = np.vstack((spine_rid[cell],np.array(spine_data[cell][idx][measure])))
            spine_rid[cell] = spine_rid[cell][1:,:]
            spine_rid[cell] = (spine_rid[cell].T/spine_rid[cell][:,0]).sum(axis=1)
            ratio_rid.append(spine_rid[cell]/dend_rid[cell])
        return np.asarray(ratio_rid)
class PlottingWidgetmRNA():
    def __init__(self,fsize=16,tsize=25,fam='Source Code Pro',pixelden = 100,lw=3.0,width=10,height=8):
        rc('font',
            family=fam,
            size=fsize)
    
        rc('figure', 
            dpi=pixelden,
            titlesize=tsize,
            titleweight='heavy',
            figsize = (width,height))
    
        rc('axes', 
            # linewidth=lw,
            titlesize=20,
            titleweight='regular',
            labelsize=fsize,
            )
    
        rc('legend',
            fontsize=fsize)
        rc('xtick',
           labelsize=fsize)
        rc('ytick',
           labelsize=fsize)
        rc('boxplot.meanprops',
           linewidth=lw,
           linestyle='--')
        rc('boxplot.boxprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.capprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.capprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.flierprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.whiskerprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.medianprops',
           linewidth=lw,
           linestyle='-')
        rc('lines',
           linewidth=lw,
           markersize=4*lw)
    def CreateFolderRecursive(self,folder):
        Path(folder).mkdir(parents=True, exist_ok=True)
    def PlotRatio(self,fractions,lab,xlab,ylab,title_string,file_name,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots()
        pos = np.arange(1,3,1)
      
        bp1 = ax.boxplot(np.transpose(fractions[0:2]),widths = 0.25,positions=pos,labels=lab[0:2],showmeans=True,meanline=True,showfliers=False,meanprops = dict(color=COLORS[0]))
        bp2 = ax.boxplot(np.transpose(fractions[2:]),widths = 0.25,positions=2+pos,labels=lab[2:],showmeans=True,meanline=True,showfliers=False,meanprops = dict(color=COLORS[1]))
        
        p_values = sp.posthoc_dunn(fractions, p_adjust = 'bonferroni')
        x_points = np.asarray((pos,2+pos)).flatten()
        pairs = np.array([[1,2],[3,4],[1,4],[2,4],[1,3],[2,3]])
        for idx,pair in enumerate(pairs):
            txt = ''
            print(p_values[pair[0]][pair[1]],pair)
            if p_values[pair[0]][pair[1]] <= 0.5:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.01:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.001:
                txt += '*'
            # breakpoint()
            y_max = 8#np.array([fractions[pair[0]-1].max(),fractions[pair[1]-1].max()]).max()
            self.AnnotateText(ax,x_points[pair[0]-1],x_points[pair[1]-1],y_max,0.1,idx*0.7,txt,'k')
        self.setBoxColors(bp1,COLORS[0])
        self.setBoxColors(bp2,COLORS[1])
        if set_axis_label == 1:
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
        # plt.ylim([0,10.2])
        # plt.title(title_string,fontsize=fsize)
        # plt.tick_params(
        #     axis='x',          # changes apply to the x-axis
        #     which='both',      # both major and minor ticks are affected
        #     bottom=False,      # ticks along the bottom edge are off
        #     top=False,         # ticks along the top edge are off
        #     labelbottom=False) # labels along the bottom edge are off
        means = []
        stds = []
        for i in range(4):
            means.append(fractions[i].mean())
            stds.append(fractions[i].std())
        # breakpoint()
        for j, line in enumerate(bp1['means']):
            x, y = line.get_xydata()[1]
            text = r'${:.2f}(\pm {:.2f})$'.format(means[j], stds[j])
            ax.annotate(text, xy=(x, y))
        for j, line in enumerate(bp2['means']):
            x, y = line.get_xydata()[1]
            text = r' ${:.2f}(\pm {:.2f})$'.format(means[2+j], stds[2+j])
            ax.annotate(text, xy=(x, y))
        fig.tight_layout()
        folder = "."
        if save_it == 1:
            self.SaveFigures(file_name)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    def AnnotateText(self,ax,x1,x2,y,h,fac,txt,color,ha='center',va='bottom'):
        print(x1,x2,y,[y+fac, y+h+fac, y+h+fac, y+fac])
        if not txt == '':
            # fac = np.abs(x2-x1)*0.5
            plt.plot([x1,x1, x2,x2], [y+fac, y+h+fac, y+h+fac, y+fac], lw=1.5, c=color)
            plt.text((x1+ x2)*0.5,y+h+fac,txt, ha=ha, va=va, color=color)
    def setBoxColors(self,bp,c):
        setp(bp['boxes'], color=c)
        setp(bp['caps'], color=c)
        setp(bp['caps'], color=c)
        setp(bp['whiskers'], color=c)
        setp(bp['whiskers'], color=c)
        setp(bp['fliers'], color=c)
        setp(bp['fliers'], color=c)
        setp(bp['medians'], color=c)
    def SaveFigures(self,filename,ext_list = [".png",".svg",".pdf"]):
        for ext in ext_list:
            plt.savefig(filename+ext,dpi=300)
if __name__ == "__main__":
    S2D_glua2 = SpineDendriteRatio("/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Rizzoli data/Analysed Kanaan/GluR2_UID-Gria2_2015-11-09/")
    dend_data_glua2,spine_data_glua2 = S2D_glua2.Readfiles(exclude_cells=[])
    S2D_glua1 = SpineDendriteRatio("/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Rizzoli data/Analysed Kanaan/GluR1_UID-Gria1_2015-11-09/")
    dend_data_glua1,spine_data_glua1 = S2D_glua1.Readfiles(exclude_cells=[])
    ratio_glua2 = S2D_glua2.GetRatios('RawIntDen', dend_data_glua2,spine_data_glua2)
    ratio_glua1 = S2D_glua1.GetRatios('RawIntDen', dend_data_glua1, spine_data_glua1)
    fractions = np.array([ratio_glua2[:,2],ratio_glua2[:,3],ratio_glua1[:,2],ratio_glua1[:,3]])
    # breakpoint()
    pw = PlottingWidgetmRNA()
    op_folder = "./helm_plots/"
    pw.CreateFolderRecursive(op_folder)
    ax_label=0;
    pw.PlotRatio(fractions, ["Confocal (N={0})".format(fractions[0].shape[0]),"STED (N={0})".format(fractions[1].shape[0]),"Confocal (N={0})".format(fractions[2].shape[0]),"STED (N={0})".format(fractions[3].shape[0])],\
                  "GluA2", "Ratio (spine/dendrite)", "Ratio of spine to dendrite fluorescent intensity", "{0}spine_to_dendrite_ratio".format(op_folder),save_it=1,set_axis_label=ax_label)
        