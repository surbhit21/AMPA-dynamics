#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 15:29:23 2021

@author: surbhitwagle
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
from PlottingWidgetAMPA import *
L = 221;    #lenght of the dendrite
class DendriteWithSpines():
    def __init__(self,D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin,eta_s,eta_c):
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
                ((self.alpha+self.Lamda_ps+self.eta_s)/self.D_s)*ps - (self.beta/self.D_s)*pc,\
                    dpc,\
                        ((self.beta+self.Lamda_pc+self.eta_c)/self.D_c - self.V_p/(self.D_c*L))*pc + (self.V_p*(1-x/L)/self.D_c)*dpc - (self.alpha/self.D_c)*ps]
        
    def bc(self,ya,yb):
        return np.array([self.D_s*ya[1] + self.Jsin, self.D_s*yb[1], self.D_c*ya[3] - self.V_p*ya[2] + self.Jcin, self.D_c*yb[3]])

    def solveModel(self):
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
        # norm_ps_dist = ps_dist/ps_dist[0]
        # norm_pc_dist = pc_dist/ps_dist[0]
        # print(len(ps_dist))
        # plt.plot(x,ps_dist,label='P_s')
        # plt.plot(x,pc_dist,label='P_c')
        # plt.show()
        self.IntegralBC(delta_x,ps_dist, pc_dist)
        return x,ps_dist,pc_dist
    
    def ParameterEffect(self,scales,D_s = None,D_c = None,V_p = None,half_life_surf = None,half_life_int = None\
                          ,alpha = None,beta = None,Jsin = None,Jcin = None,eta_s = None,eta_c = None):
        output = {};
        for ind,scale in enumerate(scales):
            current_D_s = self.D_s
            current_V_p = self.V_p 
            current_half_life_surf = self.half_life_surf
            current_Jsin = self.Jsin
            current_D_c = self.D_c
            current_half_life_int = self.half_life_int
            current_Jcin = self.Jcin
            current_alpha = self.alpha
            current_beta = self.beta
            current_eta_s = self.eta_s
            current_eta_c = self.eta_c
            # output[ind] = SingleProteinModelWithoutSpines(self.D_p,self.V_p,self.half_life,self.Jin)
            if D_s:
                current_D_s *= scale
            if V_p:
                current_V_p *= scale
            if current_half_life_surf:
                current_half_life_surf *= scale
            if Jsin:
                current_Jsin *= scale
            if D_c:
                current_D_c *= scale
            if alpha:
                current_alpha *= scale
            if current_half_life_int:
                current_half_life_int *= scale
            if Jcin:
                current_Jcin *= scale
            if beta:
                current_beta *= scale
            if eta_s:
                current_eta_s *= scale
            if eta_c:
                current_eta_c *= scale
            output[ind] = DendriteWithSpines(current_D_s, current_D_c, current_V_p, current_half_life_surf, current_half_life_int, current_alpha, \
                                             current_beta, current_Jsin, current_Jcin,current_eta_s,current_eta_c)
            # print(output[ind].__dict__)
        return output
    def IntegralBC(self,delta_x,ps,pc):
        total_pc = ((self.alpha+self.eta_s)*self.Jcin)/(self.alpha*self.eta_c + self.beta*self.eta_s + self.eta_c*self.eta_s + self.alpha* self.Lamda_pc + self.eta_s * self.Lamda_pc)
        total_ps = (self.beta*self.Jcin)/(self.alpha*self.eta_c + self.beta*self.eta_s + self.eta_c*self.eta_s + self.alpha* self.Lamda_pc + self.eta_s * self.Lamda_pc)
        # x,p= self.solveModel()
        num_total_ps = np.sum(ps)*delta_x
        num_total_pc = np.sum(pc)*delta_x;
        total_p  = total_pc + total_ps
        num_total_p = num_total_pc + num_total_ps;
        print("total surface p analytic/total surface p numeric = ",total_ps/num_total_ps)
        print("total cytoplasmic p analytic/total cytoplasmic p numeric = ",total_pc/num_total_pc)
        print("total p analytic/total p numeric = ",total_p/num_total_p)
    
    
if __name__ == '__main__':
    sim_id = "001";
    SP_model1 = DendriteWithSpines(0.45,0.05,0.1,float('inf'),4.35,0.1,0.2,0.0,0.1,0.0001,0.0006);
    x,ps_dist,pc_dist = SP_model1.solveModel()
    title_string = (r"Steady-state spatial distribution"+" \n parameters:\
       "+r" $D_s$ = %.2f, half-life-surf = %.2f, $Jsin= %.2f,  \alpha = %.2f, \eta_s$ = %.1e "+" \n"+ \
        r"$D_c = {%.2f}, V_p = {%.1e}$, half-life-int = %.2f, $Jcin= %.2f, \beta = %.2f, \eta_c$ = %.1e") \
    %( SP_model1.D_s, SP_model1.half_life_surf,SP_model1.Jsin,SP_model1.alpha,SP_model1.eta_s,\
      SP_model1.D_c,SP_model1.V_p, SP_model1.half_life_int,SP_model1.Jcin,SP_model1.beta,SP_model1.eta_c);
    lab_ps =  'Surf-AMPA-R'
    lab_pc =  'Int-AMPA-R'
    x_label = r'Dendritic distance in ($\mu$M)';
    y_label= r'Protein number';
    folder= "Figures/TwoProtein/WithUptake/";
    file_name = folder+"TwoProtein_SingleSim_{0}".format(sim_id);
    pwa = PlottingWidgetAMPA()
    pwa.CreateFolderRecursive(folder)
    pwa.PlotSingleSimTwoProtein(x, ps_dist,pc_dist, lab_ps,lab_pc, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
    title_string = (r"Steady-state spatial distribution in spines"+" \n parameters:\
       "+r" $D_s$ = %.2f, half-life-surf = %.2f, $Jsin= %.2f,  \alpha = %.2f, \eta_s$ = %.1e "+" \n"+ \
        r"$D_c = {%.2f}, V_p = {%.1e}$, half-life-int = %.2f, $Jcin= %.2f, \beta = %.2f, \eta_c$ = %.1e") \
    %( SP_model1.D_s, SP_model1.half_life_surf,SP_model1.Jsin,SP_model1.alpha,SP_model1.eta_s,\
      SP_model1.D_c,SP_model1.V_p, SP_model1.half_life_int,SP_model1.Jcin,SP_model1.beta,SP_model1.eta_c);
    lab_p_spine =  'Spine-AMPA-R'
    x_label = r'Dendritic distance in ($\mu$M)';
    y_label= r'Protein number';
    file_name = "Figures/Spine_SingleSim_TwoProtein_dist_{0}".format(sim_id);
    p_spine = SP_model1.eta_s*ps_dist + SP_model1.eta_c*pc_dist;
    pwa.PlotSingleSimSingleProtein(x, p_spine, lab_p_spine, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
    
    norm_types = ["Max","Sum"];
    for norm_type in norm_types:
        title_string = (r"Steady-state spatial distribution"+" \n parameters:\
       "+r" $D_s$ = %.2f, half-life-surf = %.2f, $Jsin= %.2f,  \alpha = %.2f, \eta_s$ = %.1e "+" \n"+ \
        r"$D_c = {%.2f}, V_p = {%.1e}$, half-life-int = %.2f, $Jcin= %.2f, \beta = %.2f, \eta_c$ = %.1e") \
        %( SP_model1.D_s, SP_model1.half_life_surf,SP_model1.Jsin,SP_model1.alpha,SP_model1.eta_s,\
          SP_model1.D_c,SP_model1.V_p, SP_model1.half_life_int,SP_model1.Jcin,SP_model1.beta,SP_model1.eta_c);
        pwa.plotNormTwo(x, ps_dist,pc_dist, lab_ps,lab_pc, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
        title_string = (r"Steady-state spatial distribution in spines"+" \n parameters:\
       "+r" $D_s$ = %.2f, half-life-surf = %.2f, $Jsin= %.2f,  \alpha = %.2f, \eta_s$ = %.1e "+" \n"+ \
        r"$D_c = {%.2f}, V_p = {%.1e}$, half-life-int = %.2f, $Jcin= %.2f, \beta = %.2f, \eta_c$ = %.1e") \
        %( SP_model1.D_s, SP_model1.half_life_surf,SP_model1.Jsin,SP_model1.alpha,SP_model1.eta_s,\
          SP_model1.D_c,SP_model1.V_p, SP_model1.half_life_int,SP_model1.Jcin,SP_model1.beta,SP_model1.eta_c);
    
        file_name = folder+"Spine_TwoProtein_SingleSim_{0}".format(sim_id);
        pwa.plotNormSingle(x, p_spine, lab_p_spine, x_label, y_label, title_string, file_name,norm_type=norm_type)
        
        scales = np.array([1/5,1/2,1,2])
        modified_objs = SP_model1.ParameterEffect(scales,V_p = 1)
        labels =  ["Vp X {%.2f}"%(scale) for scale in scales]
        title_string = "Dendrite dustribution: Effect of changing velocity of active transport";
        all_ps = np.zeros((len(scales),len(x)))
        all_pc = np.zeros((len(scales),len(x)))
        all_spine = np.zeros((len(scales),len(x)))
        for ind,key in enumerate(modified_objs.keys()):
            x,all_ps[ind],all_pc[ind] = modified_objs[key].solveModel()
            all_spine[ind] = modified_objs[key].eta_s*all_ps[ind]+modified_objs[key].eta_c*all_pc[ind]
        # print(all_op)
        file_name = folder+"TwoProtein_MultiSim_velocity_{0}".format(sim_id)
        pwa.PlotMultiSimTwoProtein(x, all_ps,all_pc,labels,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMultiTwo(x, all_ps,all_pc, labels,labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
        file_name = folder+"Spine_TwoProtein_MultiSim_velocity_{0}".format(sim_id);
        title_string = "Spine dustribution: Effect of changing velocity of active transport";
        # p_spine_multi = SP_model1.eta_p*all_op;
        pwa.PlotMultipleSim(x, all_spine, labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
        pwa.plotNormMulti(x, all_spine,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
    
    
    
    
    