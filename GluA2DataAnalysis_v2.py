#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:00:30 2022

@author: surbhitwagle
"""

from AMPA_model import *
import csv
import json
from lmfit import conf_interval, minimize,Minimizer, Parameters, Parameter, report_fit, printfuncs
import matplotlib.pyplot as plt
import numpy as np
import os
import PlotBinned 
from PlottingWidgetAMPA import *
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
from functools import reduce
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp, kruskal

import scikit_posthocs as sps
import seaborn as sns
import SNSPlottingWidget as SNSP
Lengths = np.array([25,50,75,100,150,200])
scale = 0.240
COLORS = ["#005f73","#9b2226","#CA6702","#337357"]
bin_size = scale
bins = np.arange(0, Lengths.max(), bin_size)
class GluA2DataAnalysis():
    def __init__(self,DataFolder):
        self.df = DataFolder
        
    def LoadData(self, bins,exclude_cells = []):
        """
        Assumes a folder strucutre. Give the outermost folder which contains folder for each image
       df ==> cell_i ==> Rois ==> int-GLuA2/suef-GluA2/GFP
        """
        # fig,(ax0,ax1) = plt.subplots(1, 2,figsize=(20, 12), sharey=True)
        files = os.listdir(self.df) 
        int_glua2_data = {}
        surf_glua2_data = {}
        ratio_int_surf = {}
        total_ratio_int_surf = {}
        GFP_data = {}
        soma_int_data = []
        soma_surf_data = []
        raw_int_data = {}
        raw_surf_data = {}
        for l in Lengths:
            surf_glua2_data[l] = []
            int_glua2_data[l] = []
            ratio_int_surf[l] = []
            total_ratio_int_surf[l] = []
            GFP_data[l] = []
            raw_int_data[l] = []
            raw_surf_data[l] = []
        for fdx,file in enumerate(files):
            # print(file)
            if os.path.isdir(os.path.join(self.df, file)) and file not in exclude_cells:
                int_fname = os.path.join(self.df, file+"/Rois/int-GluA2/");
                surf_fname = os.path.join(self.df, file+"/Rois/surf-GluA2/")
                GFP_fname = os.path.join(self.df, file+"/Rois/GFP/")
                int_glua2 =  os.listdir(int_fname)
                surf_glua2 =  os.listdir(surf_fname)
                # print(int_glua2,surf_glua2)
                soma_file ="soma.csv"
                measure_file = "measure.csv"
                int_glua2.remove(soma_file)
                # int_glua2.remove(measure_file)
                surf_glua2.remove(soma_file)
                # surf_glua2.remove(measure_file)
                # breakpoint()
                k = 0
                for (idx, sdx) in zip(int_glua2, surf_glua2):
                    # print(k)
                    #  reading the data from csv files and removing the 1st, title line
                    int_data =  self.ReadCSVFull(int_fname+"/"+idx)
                    surf_data =  self.ReadCSVFull(surf_fname+"/"+sdx)
                    gfp_unit_data = self.ReadCSVFull(GFP_fname+"/"+idx)
                    # breakpoint()
                    
                    # calculating the binned sum
                    binned_sum_int = self.BinnedSum(int_data,bins,0,idx)
                    binned_sum_surf = self.BinnedSum(surf_data,bins,0,sdx)
                    binned_gfp_data = self.BinnedSum(gfp_unit_data,bins,0,sdx)
                    
                    # calculating the ratio for dendrites and for soma
                    binned_sum_ratio = binned_sum_int[:,1]/binned_sum_surf[:,1]
                    total_sum_ratio = binned_sum_int[:,1].sum()/binned_sum_surf[:,1].sum()
                    
                    for l in Lengths:
                        if l <= int_data[-1,0]:
                            
                            # appending the binned ratios for appropriate lengths 
                           surf_glua2_data[l].append(binned_sum_surf[:,1]) 
                           int_glua2_data[l].append(binned_sum_int[:,1])
                           GFP_data[l].append(binned_gfp_data[:,1])
                           #appending the ratios for appropriate lengths 
                           ratio_int_surf[l].append(binned_sum_ratio) 
                           total_ratio_int_surf[l].append(total_sum_ratio)
                           
                           #getting raw values only upto length l
                           
                           x_n = int(np.ceil(l/scale))
                           raw_int_data[l].append(int_data[0:x_n,1])
                           raw_surf_data[l].append(surf_data[0:x_n,1])
                           
                soma_int_data.append(self.ReadCSVFull(int_fname+"/"+soma_file))
                soma_surf_data.append(self.ReadCSVFull(surf_fname+"/"+soma_file))
                # gfp_all_data.append(gfp_unit_data)
        for l in Lengths:
            surf_glua2_data[l] = np.asarray(surf_glua2_data[l])
            int_glua2_data[l] = np.asarray(int_glua2_data[l])
            ratio_int_surf[l] = np.asarray(ratio_int_surf[l])
            total_ratio_int_surf[l] = np.asarray(total_ratio_int_surf[l])
            raw_int_data[l] = np.asarray(raw_int_data[l])
            raw_surf_data[l] = np.asarray(raw_surf_data[l])
            GFP_data[l] = np.asanyarray(GFP_data[l])
        soma_int_data = np.asarray(soma_int_data)
        soma_surf_data = np.asarray(soma_surf_data)
        # gfp_all_data = np.asarray(gfp_all_data)
        # breakpoint()
        return int_glua2_data, surf_glua2_data,\
            ratio_int_surf,total_ratio_int_surf, \
                soma_int_data,soma_surf_data,raw_int_data,raw_surf_data,GFP_data
                        
    def ReadCSVFull(self,filename):
        csv_data = []
        with open(filename) as csvDataFile:
            # read file as csv file 
            csvReader = csv.reader(csvDataFile)
            # loop over rows
            for row in csvReader:
        
               # add cell [0] to list of dates
               csv_data.append(row)
        return np.asarray(csv_data[1:]).astype(float)
   
    def BinnedSum(self,arr,bins,num_col=-1,name= None):
        # print(name)
        if len(arr.shape) == 2:
            rr,cc = arr.shape
            binned_sum = np.zeros((len(bins),cc))
            digitized = bins.searchsorted(arr[:,num_col])
            breakpoint()
            digitized[0] = digitized[1]
            for c in range(0,cc):
                binned_sum[:,c] = np.bincount(digitized, weights=arr[:,c], minlength=len(bins))
            binned_sum[:,num_col] = bins
            return binned_sum[1:]
        else:
            print("quite not the shape",arr.shape)
            return np.zeros((len(bins),arr.shape[1]))
        
    def GetSumNormDistribution(self,data):
        sums = data.sum(axis=1)
        norm_data = np.transpose(data)/sums
        return np.transpose(norm_data)
    def GetSomaNormDistribution(self,data,index=1):
        d_vec = data[:,index]
        norm_data = data / d_vec[:,None]
        return norm_data 
    def PlotBinnedStats2P(self,x,mean_1,mean_2,sem1,sem2,num_sample,lab1,lab2,xlab,ylab,title_string,file_name,bin_size,rat=1,soma_rat=1,width=10,height=8,fsize=16,save_it = 1,fit_exp =0,set_axis_label=1):
        fig = plt.figure(figsize=(width, height))
        ax = fig.add_subplot(111) # the big subplot
        ax0 = fig.add_subplot(211)
        ax1 = fig.add_subplot(212)
        
        
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        
        if not (mean_1 == [] and sem1 == []):
            # ax0.errorbar(x+bin_size/2,mean_1,std_1,label=lab1,color=COLORS[0],marker='o',linestyle="None",markersize=0.1 )
            ax0.plot(x,mean_1,color=COLORS[0],marker=None,ls='-',label=lab1)
            CI1 = 1.96 * sem1
            ax0.fill_between(x, (mean_1-CI1), (mean_1+CI1), color=COLORS[0], alpha=.1)
            # ax0.plot(x+bin_size/2,mean_1-std_1,color=COLORS[0],marker=None,ls='--' )
        if not (mean_2 == [] and std_2 == []):
            print("in mean 2")
            # ax1.errorbar(x+bin_size/2,mean_2,std_2,label=lab2,color=COLORS[1],marker='d',linestyle="None",markersize=0.1 )
            CI2 = 1.96 * sem2
            ax1.plot(x,mean_2,color=COLORS[1],marker=None,ls='-',label=lab2 )
            ax1.fill_between(x, (mean_2-CI2), (mean_2+CI2), color=COLORS[1], alpha=.1)
            # ax1.plot(x+bin_size/2,mean_2-std_2,color=COLORS[1],marker=None,ls='--' )
        if set_axis_label ==1:
            ax.set_xlabel(xlab,fontsize=fsize)
            ax.set_ylabel(ylab,fontsize=fsize)
            fig.suptitle(title_string,fontsize=fsize)
        # plt.title(title_string,fontsize=fsize)
       
        folder = "."
        # ax0.set_ylim([0.5,1.5])
        # ax1.set_ylim([0.5,1.5])
        if fit_exp == 1:
            # breakpoint()
            # y1_fit, r1_squared = ExpFit(x,mean_1)
            # y2_fit, r2_squared = ExpFit(x,mean_2)
            # ax.plot(x+bin_size/2,y1_fit,'o--',c=COLORS[0],label=lab1+"-fit, $r^2$ =%0.2f" %(r1_squared))
            # ax.plot(x+bin_size/2,y2_fit,'o--',c=COLORS[1],label=lab2+"-fit, $r^2$ =%0.2f" %(r2_squared))
            x1,ps_mean,pc_mean,ps_rsquared,pc_rsquared,params = FitModel(x,np.stack((mean_1,mean_2)),rat,soma_rat)
            print(params)
            print(ps_rsquared,pc_rsquared)
            # x3,ps_l,pc_l,ps_l_rsquared,pc_l_rsquared,params_u = FitModel(x,np.stack((mean_1-std_1,mean_2-std_2)),rat,soma_rat,params)
            # x2,ps_u,pc_u,ps_u_rsquared,pc_u_rsquared,params_l = FitModel(x,np.stack((mean_1+std_1,mean_2+std_2)),rat,soma_rat,params)
            
            ax0.plot(x1,ps_mean,c=COLORS[0],label=lab1+"-fit")
            ax1.plot(x1,pc_mean,c=COLORS[1],label=lab2+"-fit" )
            # ax0.plot(x3,ps_l,c=COLORS[0],linestyle='--')
            # ax0.plot(x2,ps_u,c=COLORS[0],linestyle='--')
            # ax1.plot(x3,pc_l,c=COLORS[1],linestyle='--')
            # ax1.plot(x2,pc_u,c=COLORS[1],linestyle='--')
        ax0.legend(prop={'size': fsize})
        ax1.legend(prop={'size': fsize})
        # plt.legend(prop={'size': fsize})
        # breakpoint()
        # plt.show()
        fig.tight_layout()
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=300)
            plt.savefig(file_name+".svg",dpi=300)
            plt.savefig(file_name+".pdf",dpi=300)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    def PlotBinnedStats1P(self,x,mean_1,mean_2,std_1,std_2,lab1,lab2,xlab,ylab,title_string,file_name,bin_size,rat=1,soma_rat=1,width=10,height=8,fsize=16,save_it = 1,fit_exp =0,set_axis_label=1):
        fig,ax = plt.subplots(figsize=(width, height))
        
        # ax.plot(x,np.zeros(x.shape),'k--',label='=0')
        if not (mean_1 == [] and std_1 == []):
            ax.errorbar(x+bin_size/2,mean_1,std_1,label=lab1,color=COLORS[0],marker='o',linestyle="None" )
        if not (mean_2 == [] and std_2 == []):
            print("in mean 2")
            ax.errorbar(x+bin_size/2,mean_2,std_2,label=lab2,color=COLORS[1],marker='d',linestyle="None" )
        if set_axis_label ==1:
            ax.set_xlabel(xlab,fontsize=fsize)
            ax.set_ylabel(ylab,fontsize=fsize)
            # plt.title(title_string,fontsize=fsize)
            fig.suptitle(title_string,fontsize=fsize)
       
        folder = "."
        # plt.ylim([0.5,1.5])
        if fit_exp == 1:
            # breakpoint()
            # y1_fit, r1_squared = ExpFit(x,mean_1)
            # y2_fit, r2_squared = ExpFit(x,mean_2)
            # ax.plot(x+bin_size/2,y1_fit,'o--',c=COLORS[0],label=lab1+"-fit, $r^2$ =%0.2f" %(r1_squared))
            # ax.plot(x+bin_size/2,y2_fit,'o--',c=COLORS[1],label=lab2+"-fit, $r^2$ =%0.2f" %(r2_squared))
            x1,ps,pc,ps_rsquared,pc_rsquared,params = FitModel(x,np.stack((mean_1,mean_2)),rat,soma_rat)
            print(params)
            print(ps_rsquared,pc_rsquared)
            ax.plot(x1,ps,c=COLORS[0],label=lab1+"-fit, $r^2$=%0.2f" %(ps_rsquared))
            ax.plot(x1,pc,c=COLORS[1],label=lab2+"-fit, $r^2$=%0.2f" %(pc_rsquared))
        ax.legend(prop={'size': fsize})
        # ax1.legend(prop={'size': fsize})
        # plt.legend(prop={'size': fsize})
        # breakpoint()
        # plt.show()
        fig.tight_layout()
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=300)
            plt.savefig(file_name+".svg",dpi=300)
            plt.savefig(file_name+".pdf",dpi=300)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    def PlotBinnedStatsWithCI(self,x,mean_1,mean_2,std_1,std_2,num_sample,lab1,lab2,xlab,ylab,title_string,file_name,bin_size,rat=1,soma_rat=1,width=10,height=8,fsize=16,save_it = 1,fit_exp =0,set_axis_label=1):
        fig, ax = plt.subplots(figsize=(width, height))
       
        ax.plot(x,np.zeros(x.shape),'k--',label='=0')
        if not (mean_1 == [] and std_1 == []):
            # breakpoint()
            ci1 = 1.96 * std_1/np.sqrt(num_sample)
            ax.plot(x+bin_size/2,mean_1,label=lab1,color=COLORS[0],marker='o',ls='-' )
            # ax.plot(x+bin_size/2,(mean_1+ci1),color=COLORS[0],marker=None,ls="--")
            # ax.plot(x+bin_size/2,mean_1-ci1,color=COLORS[0],marker=None,ls="--")
            ax.fill_between(x+bin_size/2,(mean_1-ci1),(mean_1+ci1),color=COLORS[0],alpha=0.2)
        if not (mean_2 == [] and std_2 == []):
            ci2 = 1.96 * std_2/np.sqrt(num_sample)
            ax.plot(x+bin_size/2,mean_2,label=lab2,color=COLORS[1],marker='o',ls='-' )
            ax.fill_between(x+bin_size/2,(mean_2-ci2),(mean_2+ci2),color=COLORS[1],alpha=0.2)
            # ax.plot(x+bin_size/2,mean_2+ci2,color=COLORS[1],marker=None,ls="--")
            # ax.plot(x+bin_size/2,mean_2-ci2,color=COLORS[1],marker=None,ls="--")
        if set_axis_label ==1:
            ax.set_xlabel(xlab,fontsize=fsize)
            ax.set_ylabel(ylab,fontsize=fsize)
            # plt.title(title_string,fontsize=fsize)
            fig.suptitle(title_string,fontsize=fsize)
       
        folder = "."
       
        if fit_exp == 1:
            # breakpoint()
            # y1_fit, r1_squared = ExpFit(x,mean_1)
            # y2_fit, r2_squared = ExpFit(x,mean_2)
            # ax.plot(x+bin_size/2,y1_fit,'o--',c=COLORS[0],label=lab1+"-fit, $r^2$ =%0.2f" %(r1_squared))
            # ax.plot(x+bin_size/2,y2_fit,'o--',c=COLORS[1],label=lab2+"-fit, $r^2$ =%0.2f" %(r2_squared))
            x1,ps,pc,ps_rsquared,pc_rsquared,params = FitModel(x,np.stack((mean_1,mean_2)),rat,soma_rat)
            print(params)
            ax.plot(x1,pc,c=COLORS[0],label=lab1+"-fit, $r^2$ =%0.2f" %(pc_rsquared))
            ax.plot(x1,ps,c=COLORS[1],label=lab2+"-fit, $r^2$ =%0.2f" %(ps_rsquared))
        plt.legend(prop={'size': fsize})
        # breakpoint()
        # plt.show()
        fig.tight_layout()
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=300)
            plt.savefig(file_name+".svg",dpi=300)
            plt.savefig(file_name+".pdf",dpi=300)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
   

    def PlotRatio(self,fractions,lab,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots(figsize=(width, height))
        # plt.rc('font', **{'family':'serif','serif':['Palatino']})
        # plt.rc('text', usetex=True)
        # ax.scatter(x_data1,y_data1,label=lab1)
        # ax.scatter(x_data2,y_data2,label=lab2)
        # ax.plot(x,np.zeros(x.shape),'k--',label='=0')
        pos = np.arange(1,2,1)
        # breakpoint()
        # for i in range(0,fractions.shape[0]):
        bp1 = ax.boxplot(fractions[0],widths = 0.25,positions=pos,showmeans=True,meanline=True,labels=lab[0:1],showfliers=False)
        bp2 = ax.boxplot(fractions[1],widths = 0.25,positions=1+pos,showmeans=True,meanline=True,labels=lab[1:2],showfliers=False)
        
        p_values = sps.posthoc_dunn(fractions, p_adjust = 'bonferroni')
        x_points = np.asarray((pos,1+pos)).flatten()
        pairs = np.array([[1,2]])
        for idx,pair in enumerate(pairs):
            txt = ''
            print(p_values[pair[0]][pair[1]],pair)
            if p_values[pair[0]][pair[1]] <= 0.05:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.01:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.001:
                txt += '*'
            y_max = np.array([fractions[pair[0]-1].max(),fractions[pair[1]-1].max()]).max()
            self.AnnotateText(ax,x_points[pair[0]-1],x_points[pair[1]-1],y_max,0.01,txt,'k')
        self.setBoxColors(bp1,COLORS[0])
        self.setBoxColors(bp2,COLORS[1])
        # ax.bars(cells,dend,label=lab2,color=CB91_Blue)
        if set_axis_label ==1:
            ax.set_xlabel(xlab,fontsize=fsize)
            ax.set_ylabel(ylab,fontsize=fsize)
            # plt.title(title_string,fontsize=fsize)
            fig.suptitle(title_string,fontsize=fsize)
        # plt.legend(prop={'size': fsize})
        # plt.show()
        means = []
        stds = []
        for i in range(2):
            means.append(fractions[i].mean())
            stds.append(fractions[i].std())
        # breakpoint()
        for j, line in enumerate(bp1['means']):
            x, y = line.get_xydata()[1]
            text = r'${:.2f}(\pm {:.2f})$'.format(means[j], stds[j])
            ax.annotate(text, xy=(x, y))
        for j, line in enumerate(bp2['means']):
            x, y = line.get_xydata()[1]
            text = r' ${:.2f}(\pm {:.2f})$'.format(means[1+j], stds[1+j])
            ax.annotate(text, xy=(x, y))
        fig.tight_layout()
        folder = "."
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=300)
            plt.savefig(file_name+".svg",dpi=300)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    def AnnotateText(self,ax,x1,x2,y,h,txt,color,ha='center',va='bottom'):
        # print(x,y,txt)
        if not txt == '':
            fac = np.abs(x2-x1)*0.08
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
    
    def PlotMultiTwoProtein(self,x,ps,pc,lab_ps,lab_pc,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,adjust_x=0,adjust_y=0,save_it = 1):
        # fig, ax = plt.subplots(figsize=(width, height))
        if x.shape != (ps[0].shape) and (ps[0].shape) != (pc[0].shape):
            print("data shapes are not equal x is of size ")
        else:
            # markers = ['-','-.'];
            fig= plt.figure(figsize=(width, height))
            ax = fig.add_subplot(111) # the big subplot
            ax0 = fig.add_subplot(211)
            ax1 = fig.add_subplot(212)
            
            #  turn off the big suplot axis labels and ticks
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
            ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

            ax0.plot(x,ps,color=COLORS[0],alpha=0.1)
            ax0.plot(x,ps.mean(axis=1),color=COLORS[0],label=lab_ps)
            ax0.legend(prop={'size': fsize})
            
            ax1.plot(x,pc,color=COLORS[1],alpha=0.1)
            ax1.plot(x,pc.mean(axis=1),color=COLORS[1],label=lab_pc)
            if adjust_y == 1:
                ax0.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if adjust_x == 1:
                ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            ax.set_ylabel(ylab,fontsize=fsize)
            ax.set_xlabel(xlab,fontsize=fsize)
            
            fig.suptitle(title_string,fontsize=fsize)
            ax1.legend(prop={'size': fsize})
            # plt.show()
            folder = "."
            print(file_name)
            if save_it == 1:
                plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=300)
                plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=300)
                print("saved figures to: {%s/%s}" %(folder, file_name))
            else:
                print("Plots not saved")
            plt.show()
            
def exponential(x, a, b):
    return a*np.exp(b*x)

def GetAMPADistribution(x,a,b):
    return RunSim5(x,scale,a, b)
def ExpFit(xdata,ydata):
    param_bounds=([0,0],[np.inf,np.inf])
    popt, pcov = curve_fit(GetAMPADistribution, xdata, ydata,bounds = param_bounds)
    print(popt)
    y_fit = exponential(xdata, *popt)
    residuals = ydata- y_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return y_fit,r_squared

def FitModel(x,data,rat,soma_rat,pars=[]):
    if pars == []:
        fit_paramas = Parameters()
        np.random.seed(2022)
        dc_min = 0.05
        dc_max = 10
        dc_init = np.random.uniform(dc_min,dc_max)
        ds_min = 0.1
        ds_max = 1.0
        ds_init = np.random.uniform(ds_min,ds_max)
        vp_min = 0.0
        vp_max =1.5
        vp_init = np.random.uniform(vp_min,vp_max)
        
        # breakpoint()
        print("dc_init = ",dc_init)
        print("ds_init = ",ds_init)
        print("vp_init = ",vp_init)
        #  parameter ranges
        
        #  parameters to fit
        fit_paramas.add('dc',dc_init,min=0,max=dc_max)
        fit_paramas.add('ds',ds_init,min=0,max=ds_max)
        fit_paramas.add('vp',vp_init,min=0)
        
        # rat = 1.1 # ratio between alpha and beta
        #fixed parameters
        
        fit_paramas.add('dx',scale,vary=False)
        fit_paramas.add('half_life_surf',float('inf'),vary=False)
        fit_paramas.add('half_life_int',1.95,vary=False)
        fit_paramas.add('alpha',1,vary=False)
        fit_paramas.add('beta',1,vary=False)
        fit_paramas.add('Jsin',0.021/soma_rat,vary=False)
        fit_paramas.add('Jcin',0.021,vary=False)
        fit_paramas.add('eta_s_max',60,vary=False)
        fit_paramas.add('eta_s_zero',1/(15*60),vary=False)
        fit_paramas.add('gamma',1/(15*60),vary=False)
        fit_paramas.add('rat',rat,vary=False)
    else:
        fit_paramas = pars
    
    # breakpoint()
    # mini = Minimizer(resudual,fit_paramas,fcn_kws={'x':x, 'data':data})
    
    # out1 = minimize(resudual,fit_paramas,method='Nelder',tol=1e-10,args=(x, data))
    # report_fit(out1.params)
    # breakpoint()
    resudual(fit_paramas,x,data)
    out2 = minimize(resudual,params=fit_paramas,method='leastsq',args=(x, data))
    report_fit(out2.params)
    # breakpoint()
    # breakpoint()
    # ci, trace = conf_interval(mini, out2, sigmas=[1, 2], trace=True)
    # printfuncs.report_ci(ci)
    # breakpoint()
    
    return FittedCalculation(out2.params,x,data)
def resudual(paras,x=None,data=None):
    
    # pc_needed = pc_needed/(pc_binned.sum()*delta_x)
    # breakpoint()
    
    ps_needed,pc_needed = GetRequiredDist(paras,x,data)
    # print(x,data,ps_needed,pc_needed)
    ps_res = data[0] - ps_needed 
    pc_res = data[1] - pc_needed
    resd = np.stack((ps_res,pc_res))
    return pc_res #resd.flatten()

def GetSlidingWindowMean(data,window_len,mode='same'):
    try:
        conv_window = np.ones(window_len)*(1/window_len)
        sw_mean = np.convolve(data, conv_window,mode=mode)
        return sw_mean
    except:
        print("exception")
def GetSlidingWindowMeanMatrix(data,window_len,mode='same'):
    
    if len(data.shape) != 2:
        return ("data is not matrix ")
    print("here")
    op_matrix = []
    # op_matrix = np.ones((data.shape[0],op.shape[0]))
    
    for d in data:
        # breakpoint()
        op_matrix.append(GetSlidingWindowMean(d,window_len,mode))
    op_matrix = np.asarray(op_matrix)
    return op_matrix
def GetRequiredDist(paras,x,data):
    
    x1,ps_model, pc_model = GetParamAndModelDist(paras)
    G2DA1 = GluA2DataAnalysis("/")
    # binning the model distribution in the same size
    # ps_binned = G2DA1.BinnedSum(np.column_stack((x1,ps_model)), bins,0)[1:data.shape[1]+1,1]
    # pc_binned = G2DA1.BinnedSum(np.column_stack((x1,pc_model)), bins,0)[1:data.shape[1]+1,1]
    # taking the first N bins
    ps_needed = ps_model[0:x.shape[0]]
    # ps_needed = ps_needed/(ps_binned.sum()*delta_x)
    pc_needed = pc_model[0:x.shape[0]]
    
    # normalizing with the first bin / same as soma
    ps_needed = ps_needed/ps_needed[0]
    pc_needed = pc_needed/pc_needed[0]
    return ps_needed,pc_needed

def GetParamAndModelDist(paras):
    # reading parameters to fit
    D_c = paras['dc'].value
    D_s = paras['ds'].value
    
    
    #  reading fixed parameters
    V_p = paras['vp'].value
    delta_x = paras['dx'].value
    half_life_surf = paras['half_life_surf'].value
    half_life_int = paras['half_life_int'].value
    alpha = paras['alpha'].value
    beta = paras['beta'].value
    Jsin = paras['Jsin'].value
    Jcin = paras['Jcin'].value
    eta_s_max = paras['eta_s_max'].value
    eta_s_zero = paras['eta_s_zero'].value
    gamma = paras['gamma'].value
    rat = paras['rat'].value
    # breakpoint()
    # return model distribution
    return RunModel(D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin,eta_s_max,eta_s_zero,gamma,delta_x,rat)
    
def RunModel(D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin,eta_s_max,eta_s_zero,gamma,delta_x,rat):
    SP_model1 = DendriteWithStochasticSpinesConstantV(D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin,eta_s_max,eta_s_zero,gamma,rat);
    x1,ps_dist,pc_dist = SP_model1.solveModel(delta_x)
    # breakpoint()
    ps_sum,pc_sum = SP_model1.IntegralBC(ps_dist, pc_dist, delta_x)
    # returning sum normalized distribution
    return x1,ps_dist, pc_dist
    
def FittedCalculation(paras,x,data):
    x1,ps_dist,pc_dist = GetParamAndModelDist(paras)
    ps_needed,pc_needed = GetRequiredDist(paras,x,data)
    # GetParamAndModelDist
    delta_x = paras['dx'].value
    x_n = int(np.ceil(x[-1]/delta_x))
    # breakpoint()
    ps_rsquared = R_seq(data[0],ps_needed)
    pc_rsquared = R_seq(data[1],pc_needed)
    # breakpoint()
    return x1[0:x_n],(ps_dist/ps_dist[0])[0:x_n],(pc_dist/pc_dist[0])[0:x_n],ps_rsquared,pc_rsquared,paras
def R_seq(ydata,y_fit):
    residuals = ydata- y_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

if __name__ == "__main__":
    # G2DA = GluA2DataAnalysis("/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/single images")
    G2DA = GluA2DataAnalysis("/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/GluA2/Control/old_data")
    int_glua2_data, surf_glua2_data,ratio_int_surf,\
        total_ratio_int_surf,soma_int_data,\
            soma_surf_data, raw_int_data,raw_surf_data,GFP_data = G2DA.LoadData(bins,exclude_cells=[])
    # breakpoint()
    pwa = PlottingWidgetAMPA()
    Dir = "Data/"
    op_folder = "./"
    pwa.CreateFolderRecursive(Dir)
    win_len = 15
    # breakpoint()
    norm_int = {}
    norm_surf = {}
    mean_int = {}
    mean_surf = {}
    std_int = {}
    std_surf = {}
    sem_int = {}
    sem_surf = {}
    ratio_mean = {}
    ratio_std = {}
    ratio_sem = {}
    # f = open('Data/surf_data.json')
    # mean_surf = json.load(f)
    # f.close()
    # f = open('Data/int_data.json')
    # mean_int = json.load(f)
    # f.close()
    # x1,ps,pc,ps_rsquared,pc_rsquared,params = FitModel(np.arange(0, 100, bin_size),np.stack((mean_surf,mean_int)),rat=1.3)
    
    for dend_len in Lengths[3:4]:
        x1 = np.arange(0, dend_len, bin_size)
        offset = int(int(win_len/scale)/2)
        
        GFP_SW_data = GetSlidingWindowMeanMatrix(GFP_data[dend_len],int(win_len/scale))[:,offset:x1.shape[0]+offset]
        SURF_SW_data = GetSlidingWindowMeanMatrix(surf_glua2_data[dend_len],int(win_len/scale))[:,offset:x1.shape[0]+offset]
        INT_SW_data = GetSlidingWindowMeanMatrix(int_glua2_data[dend_len],int(win_len/scale))[:,offset:x1.shape[0]+offset]
       
        # breakpoint()
       
        GFP_SOMA_NORM_data = G2DA.GetSomaNormDistribution(GFP_SW_data,index=0)
        SURF_SOMA_NORM_data = G2DA.GetSomaNormDistribution(SURF_SW_data,index=0)
        INT_SOMA_NORM_data = G2DA.GetSomaNormDistribution(INT_SW_data,index=0)
        # breakpoint()
        GFP_SW_mean = GFP_SOMA_NORM_data.mean(axis=0)
        GFP_SW_std = GFP_SOMA_NORM_data.std(axis=0)
        GFP_SW_sem = GFP_SW_std/np.sqrt(GFP_SOMA_NORM_data.shape[0])
        
        SURF_SW_mean = SURF_SOMA_NORM_data.mean(axis=0)
        SURF_SW_std = SURF_SOMA_NORM_data.std(axis=0)
        SURF_SW_sem = SURF_SW_std/np.sqrt(SURF_SOMA_NORM_data.shape[0])
        SURF_SW_sem_n = SURF_SW_sem/np.sqrt(GFP_SW_mean)
        
        INT_SW_mean = INT_SOMA_NORM_data.mean(axis=0)
        INT_SW_std = INT_SOMA_NORM_data.std(axis=0)
        INT_SW_sem = INT_SW_std/np.sqrt(INT_SOMA_NORM_data.shape[0])
        INT_SW_sem_n = INT_SW_sem/np.sqrt(GFP_SW_mean)
        
        SURF_density = SURF_SW_mean/np.sqrt(GFP_SW_mean)
        INT_density = INT_SW_mean/GFP_SW_mean
        
        
        # ratio calculations
        
        norm_surf_data = SURF_SOMA_NORM_data/np.sqrt(GFP_SW_mean)
        norm_int_data = INT_SOMA_NORM_data/GFP_SW_mean
        
        
        
        # breakpoint()
        uni_wise_ratios = norm_int_data/norm_surf_data
        ratios_mean = uni_wise_ratios.mean(axis=0)
        ratios_std =  uni_wise_ratios.std(axis=0)
        ratios_sem = ratios_std/np.sqrt(uni_wise_ratios.shape[0])
        
        
        popt,yi_fit, ri_squared,chi_squ = SNSP.ExpFit("NormE",x1,ratios_mean,ratios_sem,0,+1,"GluA2")
        CI2 = 1.96 * ratios_sem
        fig, ax = plt.subplots()
        ax.plot(x1,ratios_mean,label='Cyto/Surf',color='b')
        ax.fill_between(x1, (ratios_mean-CI2), (ratios_mean+CI2), color='b', alpha=.1)
        plt.plot(x1,yi_fit,color='r',label="exp-fit")
        ax.set_xlabel("Dendritic distance in microns",fontsize=14)
        ax.set_ylabel("Cytoplasmic/Surface GluA2 ratio",fontsize=14)
        plt.legend()
        plt.savefig("./Figures/RatioDist.png",dpi=300)
        plt.show()
        
        # breakpoint()
        
        CI = 1.96 * GFP_SW_sem
        fig, ax = plt.subplots()
        ax.plot(x1,np.sqrt(GFP_SW_mean),label='GFP^(1/2)')
        ax.fill_between(x1, (np.sqrt(GFP_SW_mean)-CI), (np.sqrt(GFP_SW_mean)+CI), color='b', alpha=.1)
        ax.set_xlabel("Dendritic distance in microns",fontsize=14)
        ax.set_ylabel("Normalized Intensity",fontsize=14)
       
        width = 5
        x1_map2_mean = np.load("./Final_MAP2_{}_{}.npy".format(width,dend_len))
        x1_map2_sem = np.load("./Final_MAP2_sem_{}_{}.npy".format(width,dend_len))
        
        ax.plot(x1_map2_mean[0],x1_map2_mean[1],label='MAP2')
        ci1 = 1.96 * x1_map2_sem[1]
        ax.fill_between(x1_map2_mean[0], (x1_map2_mean[1]-ci1), (x1_map2_mean[1]+ci1), color='r', alpha=.1)
        plt.legend()
        plt.savefig(op_folder+"GFP_dist.png",dpi=300)
        
        plt.show()
        # np.save("../../Mo/Python-code/ampa-dynamics/Final_MAP2_std_{}_{}.npy".format(width,l1),x1_std)
        # np.save("../../Mo/Python-code/ampa-dynamics/Final_MAP2_sem_{}_{}.npy".format(width,l1),x1_sem)
        
        # # SURF_density = SURF_density
        # # INT_density = int_density
        ratios = INT_density/SURF_density
        # plt.plot(bins[1:int(dend_len/scale)],SURF_density[:int(dend_len/scale)-1],label='surf')
        # plt.plot(bins[1:int(dend_len/scale)],INT_density[:int(dend_len/scale)-1],label='int')
        # plt.plot(bins[1:int(dend_len/scale)],ratios[:int(dend_len/scale)-1],label='ratio')
        # # plt.ylim([0.,2])
        # plt.legend()
        # plt.show()
        
        breakpoint()
        # norm with GFP
        G2DA.PlotBinnedStats2P(x1,  SURF_density, INT_density,SURF_SW_sem_n,INT_SW_sem_n,SURF_SW_data.shape[0],  'surface','cytoplasmic', 'Dendritic distance (in microns)', "GFP Normalized fluoresence intensity", "GluA2 Fluorescent Distribution (N={0})".format(SURF_SW_data.shape[0]), \
                               "./Figures/GluA2_SW_Soma_norm_GFP_norm_dist_"+str(dend_len)+"_uM_without_fit", bin_size,rat = popt[0],soma_rat=ratios.mean(),save_it = 0,fit_exp=0,set_axis_label=1)
        
        
        # G2DA.PlotBinnedStats1P(x1, ratios,[], , [], 'cytoplasmic/surface', '', 'Dendritic distance (in microns)', "cytoplasmic/surface ratio", "Distribution of cytoplasmic/surface GluA2 fluorescent intensity (N={0})".format(ratio_int_surf[l].shape[0]), \
        #                       op_folder+"GluA2_cyto_surf_ratio_dist_"+str(l)+"_uM", bin_size,save_it = 1,fit_exp=0,set_axis_label=ax_label)
           
        # norm with MAP2
        # min_size = min(x1_map2_mean.shape[1],SURF_SW_mean.shape[0])
        # SURF_density_MAP2 = SURF_SW_mean[0:min_size]/x1_map2_mean[1,0:min_size]
        # INT_density_MAP2  = INT_SW_mean[0:min_size]/x1_map2_mean[1,0:min_size]
        # SURF_SW_sem_MAp2 = SURF_SW_sem[0:min_size]/x1_map2_mean[1,0:min_size]
        # INT_SW_sem_MAP2 = INT_SW_sem[0:min_size]/x1_map2_mean[1,0:min_size]
        
        # G2DA.PlotBinnedStats2P(x1,  SURF_density_MAP2, INT_density_MAP2,SURF_SW_sem_MAp2,INT_SW_sem_MAP2,SURF_SW_data.shape[0],  'surface','cytoplasmic', 'Dendritic distance (in microns)', "MAP2 Normalized fluoresence intensity", "GluA2 Fluorescent Distribution (N={0})".format(SURF_SW_data.shape[0]), \
        #                        op_folder+"GluA2_SW_Soma_norm_MAP2_norm_dist_"+str(dend_len)+"_uM_with_fit", bin_size,rat = ratios.mean(),soma_rat=ratios.mean(),save_it = 1,fit_exp=1,set_axis_label=1)
 
 