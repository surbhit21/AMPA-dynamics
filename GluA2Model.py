#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 14:50:04 2021

@author: surbhitwagle
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
from PlottingWidgetAMPA import *
L = 221;    #lenght of the dendrite
class GluA2Model():
    def __init__(self,D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin,eta_s,eta_c,k_c):
        self.D_s = D_s   # in uM^2/s
        self.D_c = D_c   # in uM^2/s
        self.V_p = V_p    # in uM/s
        self.half_life_surf = half_life_surf # in days
        self.Lamda_ps = np.log(2)/(self.half_life_surf*24*60*60);
        self.half_life_int = half_life_int # in days
        self.Lamda_pc = np.log(2)/(self.half_life_int*24*60*60);
        self.alpha = alpha;
        self.beta = beta;
        self.Jsin = Jsin;
        self.Jcin = Jcin;
        self.eta_s = eta_s;
        self.eta_c = eta_c;
        self.k_c = k_c;
        
    def updateModelParams(self,D_s = None,D_c = None,V_p = None,half_life_surf = None,half_life_int = None\
                          ,alpha = None,beta = None,Jsin = None,Jcin = None,eta_s = None,eta_c = None):
        if D_s:
            self.D_s = D_s   # in uM^2/s
        if D_c:
            self.D_c = D_c   # in uM^2/s
        if V_p:
            self.V_p = V_p    # in uM/s
        if half_life_surf:
            self.half_life_surf = half_life_surf # in days
            self.Lamda_ps = np.log(2)/(self.half_life_surf*24*60*60);
        if half_life_int:
            self.half_life_int = half_life_int # in days
            self.Lamda_pc = np.log(2)/(self.half_life_int*24*60*60);
        if alpha:
            self.alpha = alpha;
        if beta:
            self.beta = beta;
        if Jsin:
            self.Jsin = Jsin;
        if Jcin:
            self.Jcin = Jcin;
        if eta_s:
            self.eta_s = eta_s;
        if eta_c:
            self.eta_c = eta_c;
    
    def fun(self,x,y):
        ps,dps,pc,dpc = y
        return [dps,\
                ((self.alpha+self.Lamda_ps+self.eta_s)/self.D_s)*ps - (self.beta*np.exp(-self.k_c*x)/self.D_s)*pc,\
                    dpc,\
                        ((self.beta*np.exp(-self.k_c*x)+self.Lamda_pc+self.eta_c)/self.D_c - self.V_p/(self.D_c*L))*pc + (self.V_p*(1-x/L)/self.D_c)*dpc - (self.alpha/self.D_c)*ps]
        
    def bc(self,ya,yb):
        return np.array([self.D_s*ya[1] + self.Jsin, self.D_s*yb[1], self.D_c*ya[3] - self.V_p*ya[2] + self.Jcin, self.D_c*yb[3]])

    def solveModel(self,sim_id):
        delta_x = 0.0114; #step size for the numerical simulation
        # sim_id="001";
        # print(len(x_line))
        # solving model
        # D_p * p'' - (V_p(x)*p)' - Lamda_p*p = 0
        x=np.arange(0,L,delta_x)
        print(x)
        params=np.array([L]);
        y = np.zeros((4,x.size))
        soln = solve_bvp(self.fun, self.bc, x, y)
        print(len(soln))
        ps_dist = soln.sol(x)[0]
        pc_dist = soln.sol(x)[2]
        norm_ps_dist = ps_dist/ps_dist[0]
        norm_pc_dist = pc_dist/ps_dist[0]
        print(len(ps_dist))
        # plt.plot(x,ps_dist,label='P_s')
        # plt.plot(x,pc_dist,label='P_c')
        # plt.show()
        title_string = (r"Steady-state spatial distribution"+" \n parameters:\
           "+r" $D_s$ = %.2f, half-life-surf = %.2f, $Jsin= %.2f,  \alpha = %.2f, \eta_s$ = %.1e "+" \n"+ \
            r"$D_c = {%.2f}, V_p = {%.1e}$, half-life-int = %.2f, $Jcin= %.2f, \beta = %.2f, \eta_c$ = %.1e") \
        %( self.D_s, self.half_life_surf,self.Jsin,self.alpha,self.eta_s,\
          self.D_c,self.V_p, self.half_life_int,self.Jcin,self.beta,self.eta_c);
        lab_ps =  'Surf-AMPA-R'
        lab_pc =  'Int-AMPA-R'
        x_label = r'Dendritic distance in ($\mu$M)';
        y_label= r'Protein number';
        file_name = "TwoProtein_SingleSim_withSpines_{0}".format(sim_id);
        pwa = PlottingWidgetAMPA()
        pwa.PlotSingleSimTwoProtein(x, ps_dist,pc_dist, lab_ps,lab_pc, x_label, y_label, title_string, file_name,fsize=14,save_it = 0)

if __name__ == '__main__':
    GluA2_model1 = GluA2Model(0.45,0.05,0.1,float('inf'),4.35,0.1,0.2,0.0,0.1,0.0001,0.0006,0.001);
    GluA2_model1.solveModel("001")