#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 14:01:16 2021

@author: surbhitwagle
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
# from matplotlib import rc
from PlottingWidgetAMPA import *
L = 221     #length of the dendrite

class SingleProteinModelWihtCappedSpinesUptake():
    def __init__(self,D_p,V_p,half_life,Jin,eta_p_max,eta_p_zero):
        self.D_p = D_p   # in uM^2/s
        self.V_p = V_p    # in uM/s
        self.half_life =half_life # in days
        self.Lamda_p = np.log(2)/(self.half_life*24*60*60);
        self.Jin = Jin;
        self.eta_p_max = eta_p_max;
        self.eta_p_zero = eta_p_zero
        
    def updateModelParams(self,D_p=None,V_p=None,half_life = None,Jin = None,eta_p_max=None,eta_p_zero=None):
        if D_p:
            self.D_p = D_p
        if V_p:
            self.V_p = V_p
        if half_life:
            self.half_life = half_life
            self.Lamda_p = np.log(2)/(self.half_life*24*60*60);
        if Jin:
            self.Jin = Jin
        if eta_p_max:
            self.eta_p_max = eta_p_max
        if eta_p_zero:
            self.eta_p_zero = eta_p_zero
    
    def fun(self,x,y):
        # print(L)
        return y[1], ((self.Lamda_p- self.V_p/L)/self.D_p)*y[0] + self.eta_p_max*np.tanh(self.eta_p_zero*y[0])+ ((self.V_p*(1-x/L))/self.D_p)*y[1]
        
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
        title_string = "Steady-state spatial distribution \n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f, $\eta_{pmax}$ =%.1e,$\eta_{p0}$ =%.1e" \
        %( self.D_p, self.V_p, self.half_life,self.Jin,self.eta_p_max,self.eta_p_zero);
        lab =  'AMPA-R'
        x_label = r'Dendritic distance in ($\mu$M)';
        y_label= r'Protein number';
        file_name = "Figures/SingleProtein_SingleSim__CappedSpineUptake{0}".format(sim_id);
        pwa = PlottingWidgetAMPA()
        pwa.PlotSingleSimSingleProtein(x, p_dist, lab, x_label, y_label, title_string, file_name)
        
        title_string = "Steady-state spatial distribution \n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f, $\eta_{pmax}$ =%.1e,$\eta_{p0}$ =%.1e" \
        %( self.D_p, self.V_p, self.half_life,self.Jin,self.eta_p_max,self.eta_p_zero);
        lab_p_spine =  'Spine-AMPA-R'
        x_label = r'Dendritic distance in ($\mu$M)';
        y_label= r'Protein number';
        file_name = "Figures/Spine_SingleSim_OneProtein_CappedUptake_dist_{0}".format(sim_id);
        p_spine = self.eta_p_max*np.tanh(self.eta_p_zero*p_dist);
        pwa.PlotSingleSimSingleProtein(x, p_spine, lab_p_spine, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
        
    
if __name__ == '__main__':
    SP_model1 = SingleProteinModelWihtCappedSpinesUptake(0.45,0.1,4.35,1.0,6e-4,0.1);
    SP_model1.solveModel("001")
    # SP_model1.updateModelParams(V_p=0.001);
    # SP_model1.solveModel()