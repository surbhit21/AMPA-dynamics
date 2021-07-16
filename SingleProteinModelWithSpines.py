#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 14:25:19 2021

@author: surbhitwagle
"""
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
# from matplotlib import rc

L = 221     #length of the dendrite

class SingleProteinModelWihtSpines():
    def __init__(self,D_p,V_p,half_life,Jin,eta_p):
        self.D_p = D_p   # in uM^2/s
        self.V_p = V_p    # in uM/s
        self.half_life =half_life # in days
        self.Lamda_p = np.log(2)/(self.half_life*24*60*60);
        self.Jin = Jin;
        self.eta_p = eta_p;
        
    def updateModelParams(self,D_p=None,V_p=None,half_life = None,Jin = None,eta_p=None):
        if D_p:
            self.D_p = D_p
        if V_p:
            self.V_p = V_p
        if half_life:
            self.half_life = half_life
            self.Lamda_p = np.log(2)/(self.half_life*24*60*60);
        if Jin:
            self.Jin = Jin
        if eta_p:
            self.eta_p = eta_p
    
    def fun(self,x,y):
        # print(L)
        return y[1], ((self.Lamda_p+self.eta_p- self.V_p/L)/self.D_p)*y[0] + ((self.V_p*(1-x/L))/self.D_p)*y[1]
        
    def bc(self,ya,yb):
        return np.array([self.D_p*ya[1] - self.V_p*ya[0] + self.Jin, self.D_p*yb[1]])

    def solveModel(self,sim_id):
        delta_x = 0.0114; #step size for the numerical simulation
        # print(len(x_line))
        # solving model
        # D_p * p'' - (V_p(x)*p)' - Lamda_p*p = 0
        x=np.arange(0,L,delta_x)
        # print(x)
        # params=np.array([L]);
        y = np.zeros((2,x.size))
        soln = solve_bvp(self.fun, self.bc, x, y)
        # print(len(soln))
        p_dist = soln.sol(x)[0]
        # norm_p_dist = p_dist/p_dist[0]
        # print(len(p_dist))
        title_string = "Steady-state spatial distribution \n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f, $\eta_p$ =%.1e" \
        %( self.D_p, self.V_p, self.half_life,self.Jin,self.eta_p);
        lab =  'AMPA-R'
        x_label = r'Dendritic distance in ($\mu$M)';
        y_label= r'Protein number';
        file_name = "SingleProtein_MultiSim_{0}".format(sim_id);
        self.PlotSingleSim(x, p_dist, lab, x_label, y_label, title_string, file_name)
        # self.PlotMultipleSim( y, y, np.array([lab,lab]), x_label, y_label, title_string, file_name)
    # def SolveMultiModel(self,D_p_arr)
    def PlotSingleSim(self,x,y,lab,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
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
    
    
    def PlotMultipleSim(self,x,y,lab,xlab,ylab,title_string,file_name,width=8,height=8,fsize=16,save_it = 1):
        # plotting multiple simulations
        # useful when trying single parameter manipulations
        
        if x.shape != y.shape:
            print("data shapes are not equal x is of size ", x.shape , "and y is of size ",y.shape )
        else:
            fig, ax = plt.subplots(figsize=(width, height))
            plt.rc('font', **{'family':'serif','serif':['Palatino']})
            plt.rc('text', usetex=True)
            y_shape = y.shape
            print(y_shape[0])
            for i in range(0,y_shape[0]):
                ax.plot(x[i],y[i],label=lab[i])
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
        
    
if __name__ == '__main__':
    SP_model1 = SingleProteinModelWihtSpines(0.45,0.1,4.35,0.1,6e-4);
    SP_model1.solveModel("001")
    # SP_model1.updateModelParams(V_p=0.001);
    # SP_model1.solveModel()