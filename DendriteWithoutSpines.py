#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 11:52:04 2021

@author: surbhitwagle
"""
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
from PlottingWidgetAMPA import *
L = 221;    #lenght of the dendrite
class DendriteWithoutSpines():
    def __init__(self,D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin):
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
        
    def updateModelParams(self,D_s = None,D_c = None,V_p = None,half_life_surf = None,half_life_int = None\
                          ,alpha = None,beta = None,Jsin = None,Jcin = None):
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
    
    def fun(self,x,y):
        ps,dps,pc,dpc = y
        return [dps,\
                ((self.alpha+self.Lamda_ps)/self.D_s)*ps - (self.beta/self.D_s)*pc,\
                    dpc,\
                        ((self.beta+self.Lamda_pc)/self.D_c - self.V_p/(self.D_c*L))*pc + (self.V_p*(1-x/L)/self.D_c)*dpc - (self.alpha/self.D_c)*ps]
        
    def bc(self,ya,yb):
        return np.array([self.D_s*ya[1] + self.Jsin, self.D_s*yb[1], self.D_c*ya[3] - self.V_p*ya[2] + self.Jcin, self.D_c*yb[3]])

    def solveModel(self):
        delta_x = 0.0114; #step size for the numerical simulation
        # solving model
        # D_p * p'' - (V_p(x)*p)' - Lamda_p*p = 0
        x=np.arange(0,L,delta_x)
        # print(x)
        params=np.array([L]);
        y = np.zeros((4,x.size))
        soln = solve_bvp(self.fun, self.bc, x, y)
        # print(len(soln))
        ps_dist = soln.sol(x)[0]
        pc_dist = soln.sol(x)[2]
        # print(len(ps_dist))
        return x,ps_dist,pc_dist
        
    
    def ParameterEffect(self,scales,D_s = None,D_c = None,V_p = None,half_life_surf = None,half_life_int = None\
                          ,alpha = None,beta = None,Jsin = None,Jcin = None):
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
            output[ind] = DendriteWithoutSpines(current_D_s, current_D_c, current_V_p, current_half_life_surf, current_half_life_int, current_alpha, current_beta, current_Jsin, current_Jcin)
            # print(output[ind].__dict__)
        return output

    def IntegralBC(self):
        total_p = (self.Jcin/self.Lamda_pc)*(1+self.beta/self.alpha);# + self.Jsin/self.Lamda_ps
        x,ps,pc = self.solveModel()
        num_total_p = (np.sum(ps) + np.sum(pc))*(x[1]-x[0])
        print("total p analytic/total p numeric = ",total_p/num_total_p)
if __name__ == '__main__':
    sim_id = "001";
    SP_model1 = DendriteWithoutSpines(0.45,0.05,0.1,float('inf'),4.35,0.1,0.2,0.0,0.1);
    x,ps_dist,pc_dist=SP_model1.solveModel()
    folder= "Figures/TwoProtein/NoUptake/";
    file_name = folder+"SingleProtein_SingleSim_{0}".format(sim_id);
    pwa = PlottingWidgetAMPA()
    pwa.CreateFolderRecursive(folder)
    title_string = (r"Steady-state spatial distribution"+" \n parameters:\
           "+r" $D_s$ = %.2f, half-life-surf = %.2f, $Jsin= %.2f,  \alpha $ = %.2f"+" \n"+ \
            r"$D_c = {%.2f}, V_p = {%.1e}$, half-life-int = %.2f, $Jcin= %.2f, \beta$ = %.2f") \
        %( SP_model1.D_s, SP_model1.half_life_surf,SP_model1.Jsin,SP_model1.alpha, SP_model1.D_c,SP_model1.V_p, SP_model1.half_life_int,SP_model1.Jcin,SP_model1.beta);
    lab_ps =  'Surf-AMPA-R'
    lab_pc =  'Int-AMPA-R'
    x_label = r'Dendritic distance in ($\mu$M)';
    y_label= r'Protein number';
    file_name = folder+"TwoProtein_SingleSim_withoutSpines_{0}".format(sim_id);
    pwa = PlottingWidgetAMPA()
    norm_type = "Max"
    pwa.PlotSingleSimTwoProtein(x, ps_dist,pc_dist, lab_ps,lab_pc, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
    pwa.plotNormTwo(x, ps_dist,pc_dist, lab_ps,lab_pc, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
    SP_model1.IntegralBC()
    scales = np.array([1/5,1/2,1,2])
    modified_objs = SP_model1.ParameterEffect(scales,V_p = 1)
    labels =  ["Vp X {%.2f}"%(scale) for scale in scales]
    title_string = "Effect of changing velocity of active transport";
    all_ps = np.zeros((len(scales),len(x)))
    all_pc = np.zeros((len(scales),len(x)))
    for ind,key in enumerate(modified_objs.keys()):
        _,all_ps[ind],all_pc[ind] = modified_objs[key].solveModel()
        modified_objs[key].IntegralBC()
    # print(all_op)
    file_name = folder+"TwoProtein_MultiSim_velocity_{0}".format(sim_id)
    pwa.PlotMultiSimTwoProtein(x, all_ps,all_pc,labels,labels,x_label,y_label,title_string,file_name,save_it=1)
    pwa.plotNormMultiTwo(x, all_ps,all_pc, labels,labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
    
    scales = np.array([1/100,1,100,1000])
    title_string = "Effect of changing surface diffusion constant";
    modified_objs = SP_model1.ParameterEffect(scales,D_s = 1)
    labels =  ["Dp X {%.2f}"%(scale) for scale in scales]
    # all_op = np.zeros((len(scales),len(x)))
    for ind,key in enumerate(modified_objs.keys()):
        _,all_ps[ind],all_pc[ind] = modified_objs[key].solveModel()
        modified_objs[key].IntegralBC()
    # print(all_op)
    file_name = folder+"TwoProtein_MultiSim_surface_Diffusion_{0}".format(sim_id)
    pwa.PlotMultiSimTwoProtein(x, all_ps,all_pc,labels,labels,x_label,y_label,title_string,file_name,save_it=1)
    pwa.plotNormMultiTwo(x, all_ps,all_pc, labels,labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
    
    scales = np.array([1/5,1/2,1,2])
    title_string = "Effect of changing half-life-int";
    modified_objs = SP_model1.ParameterEffect(scales,half_life_int = 1)
    labels =  ["T-half-int X {%.2f}"%(scale) for scale in scales]
    # all_op = np.zeros((len(scales),len(x)))
    for ind,key in enumerate(modified_objs.keys()):
        _,all_ps[ind],all_pc[ind] = modified_objs[key].solveModel()
        modified_objs[key].IntegralBC()
    # print(all_op)
    file_name = folder+"TwoProtein_MultiSim_half_life_{0}".format(sim_id)
    pwa.PlotMultiSimTwoProtein(x, all_ps,all_pc,labels,labels,x_label,y_label,title_string,file_name,save_it=1)
    pwa.plotNormMultiTwo(x, all_ps,all_pc, labels,labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
    
    scales = np.array([1/5,1/2,1,2])
    title_string = "Effect of changing endocytosis/exocytosis rates";
    modified_objs = SP_model1.ParameterEffect(scales,alpha = 1,beta = 1)
    labels =  ["alpha,beta X {%.2f}"%(scale) for scale in scales]
    # all_op = np.zeros((len(scales),len(x)))
    for ind,key in enumerate(modified_objs.keys()):
        _,all_ps[ind],all_pc[ind] = modified_objs[key].solveModel()
        modified_objs[key].IntegralBC()
    # print(all_op)
    file_name = folder+"TwoProtein_MultiSim_recycling_{0}".format(sim_id)
    pwa.PlotMultiSimTwoProtein(x, all_ps,all_pc,labels,labels,x_label,y_label,title_string,file_name,save_it=1)
    pwa.plotNormMultiTwo(x, all_ps,all_pc, labels,labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
    
    scales = np.array([1/5,1/2,1,2])
    title_string = "Effect of changing endocytosis rates";
    modified_objs = SP_model1.ParameterEffect(scales,alpha = 1)
    labels =  ["alpha X {%.2f}"%(scale) for scale in scales]
    # all_op = np.zeros((len(scales),len(x)))
    for ind,key in enumerate(modified_objs.keys()):
        _,all_ps[ind],all_pc[ind] = modified_objs[key].solveModel()
    # print(all_op)
    file_name = folder+"TwoProtein_MultiSim_endocytosis_{0}".format(sim_id)
    pwa.PlotMultiSimTwoProtein(x, all_ps,all_pc,labels,labels,x_label,y_label,title_string,file_name,save_it=1)
    pwa.plotNormMultiTwo(x, all_ps,all_pc, labels,labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
    
    scales = np.array([1/5,1/2,1,2])
    title_string = "Effect of changing exocytosis rates";
    modified_objs = SP_model1.ParameterEffect(scales,beta = 1)
    labels =  ["beta X {%.2f}"%(scale) for scale in scales]
    # all_op = np.zeros((len(scales),len(x)))
    for ind,key in enumerate(modified_objs.keys()):
        _,all_ps[ind],all_pc[ind] = modified_objs[key].solveModel()
        modified_objs[key].IntegralBC()
    # print(all_op)
    file_name = folder+"TwoProtein_MultiSim_exocytosis_{0}".format(sim_id)
    pwa.PlotMultiSimTwoProtein(x, all_ps,all_pc,labels,labels,x_label,y_label,title_string,file_name,save_it=1)
    pwa.plotNormMultiTwo(x, all_ps,all_pc, labels,labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1,norm_type=norm_type)
    