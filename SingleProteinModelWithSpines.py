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
from PlottingWidgetAMPA import *
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

    def solveModel(self):
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
        self.IntegralBC(delta_x, p_dist)
        return x,p_dist

    def ParameterEffect(self,scales,D_p=None,V_p = None,half_life = None,Jin = None,eta_p=None):
        output = {};
        for ind,scale in enumerate(scales):
            current_D_p = self.D_p
            current_V_p = self.V_p 
            current_half_life = self.half_life
            current_Jin = self.Jin
            current_eta_p = self.eta_p
            # output[ind] = SingleProteinModelWithoutSpines(self.D_p,self.V_p,self.half_life,self.Jin)
            if D_p:
                current_D_p *= scale
            if V_p:
                current_V_p *= scale
            if half_life:
                current_half_life *= scale
            if Jin:
                current_Jin *= scale
            if eta_p:
                current_eta_p *= scale
            output[ind] = SingleProteinModelWihtSpines(current_D_p,current_V_p,current_half_life,current_Jin,current_eta_p)
            print(output[ind].__dict__)
        return output

    def IntegralBC(self,delta_x,p):
        total_p = self.Jin/(self.Lamda_p+self.eta_p);
        # x,p= self.solveModel()
        num_total_p = np.sum(p)*delta_x
        print("total p analytic/total p numeric = ",total_p/num_total_p)
    
if __name__ == '__main__':
    sim_id = "001"
    SP_model1 = SingleProteinModelWihtSpines(0.45,0.1,4.35,1.0,6e-4);
    x,p_dist = SP_model1.solveModel()
    # SP_model1.IntegralBC(x,p_dist)
    # SP_model1.updateModelParams(V_p=0.001);
    # SP_model1.solveModel()
    title_string = "Steady-state spatial distribution \n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f, $\eta_p$ =%.1e" \
        %( SP_model1.D_p, SP_model1.V_p, SP_model1.half_life,SP_model1.Jin,SP_model1.eta_p);
    lab =  'AMPA-R'
    x_label = r'Dendritic distance (in $\mu$M)';
    y_label= r'Protein number';
    folder= "Figures/OneProtein/WithUptake/";
    file_name = folder+"SingleProtein_SingleSim_{0}".format(sim_id);
    
    pwa = PlottingWidgetAMPA()
    pwa.CreateFolderRecursive(folder)
    pwa.PlotSingleSimSingleProtein(x, p_dist, lab, x_label, y_label, title_string, file_name)
    
    title_string = "Steady-state spatial distribution in Spines\n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f, $\eta_p$ =%.1e" \
    %( SP_model1.D_p, SP_model1.V_p, SP_model1.half_life,SP_model1.Jin,SP_model1.eta_p);
    lab_p_spine =  'Spine-AMPA-R'
    x_label = r'Dendritic distance (in $\mu$M)';
    y_label= r'Protein number';
    folder= "Figures/OneProtein/WithUptake/";
    file_name = folder+"Spine_SingleProtein_SingleSim_{0}".format(sim_id);
    p_spine = SP_model1.eta_p*p_dist;
    
    pwa.PlotSingleSimSingleProtein(x, p_spine, lab_p_spine, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
    
    
    norm_types = ["Max","Sum"];
    for norm_type in norm_types:
        title_string = "Steady-state spatial distribution \n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f, $\eta_p$ =%.1e" \
            %( SP_model1.D_p, SP_model1.V_p, SP_model1.half_life,SP_model1.Jin,SP_model1.eta_p);
        file_name = folder+"SingleProtein_SingleSim_{0}".format(sim_id);
        pwa.plotNormSingle(x, p_dist, lab, x_label, y_label, title_string, file_name)
        title_string = "Steady-state spatial distribution in Spines\n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f, $\eta_p$ =%.1e" \
            %( SP_model1.D_p, SP_model1.V_p, SP_model1.half_life,SP_model1.Jin,SP_model1.eta_p);
        file_name = folder+"Spine_SingleProtein_SingleSim_{0}".format(sim_id);
        pwa.plotNormSingle(x, p_spine, lab, x_label, y_label, title_string, file_name,norm_type=norm_type)
        
        scales = np.array([1/5,1/2,1,2])
        modified_objs = SP_model1.ParameterEffect(scales,V_p = 1)
        labels =  ["Vp X {%.2f}"%(scale) for scale in scales]
        title_string = "Dendrite dustribution: Effect of changing velocity of active transport";
        
        for ind,key in enumerate(modified_objs.keys()):
            x,all_op[ind] = modified_objs[key].solveModel()
            # modified_objs[key].IntegralBC(x,all_op[ind])
            all_sp[ind] = modified_objs[key].eta_p*all_op[ind]
        # print(all_op)
        file_name = folder+"SingleProtein_MultiSim_velocity_{0}".format(sim_id)
        pwa.PlotMultipleSim(x, all_op,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMulti(x, all_op,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        file_name = folder+"Spine_SingleProtein_MultiSim_velocity_{0}".format(sim_id);
        title_string = "Spine dustribution: Effect of changing velocity of active transport";
        # p_spine_multi = SP_model1.eta_p*all_op;
        pwa.PlotMultipleSim(x, all_sp, labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
        pwa.plotNormMulti(x, all_sp,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        
        scales = np.array([1/100,1,100,1000])
        title_string = "Dendrite dustribution: Effect of changing diffusion constant";
        modified_objs = SP_model1.ParameterEffect(scales,D_p = 1)
        labels =  ["Dp X {%.2f}"%(scale) for scale in scales]
        all_op = np.zeros((len(scales),len(x)))
        all_sp = np.zeros((len(scales),len(x)))
        for ind,key in enumerate(modified_objs.keys()):
            x,all_op[ind] = modified_objs[key].solveModel()
            all_sp[ind] = modified_objs[key].eta_p*all_op[ind]
            # modified_objs[key].IntegralBC(x,all_op[ind])
        print(all_op)
        file_name = folder+"SingleProtein_MultiSim_Diffusion_{0}".format(sim_id)
        pwa.PlotMultipleSim(x, all_op,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMulti(x, all_op,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        file_name = folder+"Spine_SingleProtein_MultiSim_Diffusion_{0}".format(sim_id);
        title_string = "Spine dustribution: Effect of changing velocity of active transport";
        # p_spine_multi = SP_model1.eta_p*all_op;
        pwa.PlotMultipleSim(x, all_sp, labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
        pwa.plotNormMulti(x, all_sp,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        
        scales = np.array([1/5,1/2,1,2])
        title_string = "Dendrite dustribution: Effect of changing half-life";
        modified_objs = SP_model1.ParameterEffect(scales,half_life = 1)
        labels =  ["T-half X {%.2f}"%(scale) for scale in scales]
        all_op = np.zeros((len(scales),len(x)))
        all_sp = np.zeros((len(scales),len(x)))
        for ind,key in enumerate(modified_objs.keys()):
            x,all_op[ind] = modified_objs[key].solveModel()
            all_sp[ind] = modified_objs[key].eta_p*all_op[ind]
            # modified_objs[key].IntegralBC(x,all_op[ind])
        print(all_op)
        file_name = folder+"SingleProtein_MultiSim_half_life_{0}".format(sim_id)
        pwa.PlotMultipleSim(x, all_op,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMulti(x, all_op,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        file_name = folder+"Spine_SingleProtein_MultiSim_half_life_{0}".format(sim_id);
        title_string = "Spine dustribution: Effect of changing half-life of protein";
        # p_spine_multi = SP_model1.eta_p*all_op;
        pwa.PlotMultipleSim(x, all_sp, labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
        pwa.plotNormMulti(x, all_sp,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        
        scales = np.array([1/5,1/2,1,2])
        title_string = "Dendrite dustribution: Effect of changing Uptake rate";
        modified_objs = SP_model1.ParameterEffect(scales,eta_p = 1)
        labels =  ["etap X {%.2f}"%(scale) for scale in scales]
        all_op = np.zeros((len(scales),len(x)))
        all_sp = np.zeros((len(scales),len(x)))
        for ind,key in enumerate(modified_objs.keys()):
            x,all_op[ind] = modified_objs[key].solveModel()
            all_sp[ind] = modified_objs[key].eta_p*all_op[ind]
            # modified_objs[key].IntegralBC(x,all_op[ind])
        print(all_op)
        file_name = folder+"SingleProtein_MultiSim_synaptic_uptake_{0}".format(sim_id)
        pwa.PlotMultipleSim(x, all_op,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMulti(x, all_op,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        file_name = folder+"Spine_SingleProtein_MultiSim_synaptic_uptake__{0}".format(sim_id);
        title_string = "Spine dustribution: Effect of changing Uptake rate";
        # p_spine_multi = SP_model1.eta_p*all_op;
        pwa.PlotMultipleSim(x, all_sp, labels, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
        pwa.plotNormMulti(x, all_sp,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)