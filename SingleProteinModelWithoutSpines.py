#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:57:48 2021

@author: surbhitwagle
python package to solve the protein dynamics model
"""
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
# from matplotlib import rc
from  PlottingWidgetAMPA import *
L = 221     #length of the dendrite

class SingleProteinModelWithoutSpines():
    def __init__(self,D_p,V_p,half_life,Jin):
        self.D_p = D_p   # in uM^2/s
        self.V_p = V_p    # in uM/s
        self.half_life =half_life # in days
        self.Lamda_p = np.log(2)/(self.half_life*24*60*60);
        self.Jin = Jin;
        
    def updateModelParams(self,D_p=None,V_p=None,half_life = None,Jin = None):
        if D_p:
            self.D_p = D_p
        if V_p:
            self.V_p = V_p
        if half_life:
            self.half_life = half_life
            self.Lamda_p = np.log(2)/(self.half_life*24*60*60);
        if Jin:
            self.Jin = Jin
    
    def fun(self,x,y):
        # print(L)
        return y[1], ((self.Lamda_p - self.V_p/L)/self.D_p)*y[0] + ((self.V_p*(1-x/L))/self.D_p)*y[1]
        
    def bc(self,ya,yb):
        return np.array([self.D_p*ya[1] - self.V_p*ya[0] + self.Jin, self.D_p*yb[1]])

    def solveModel(self):
        delta_x = 0.0114; #step size for the numerical simulation
        # print(len(x_line))
        # solving model
        # D_p * p'' - (V_p(x)*p)' - Lamda_p*p = 0
        x=np.arange(0,L,delta_x)
        
        y = np.zeros((2,x.size))
        soln = solve_bvp(self.fun, self.bc, x, y)
        # print(len(soln))
        p_dist = soln.sol(x)[0]
        # norm_p_dist = p_dist/p_dist[0]
        # print(len(p_dist))
        self.IntegralBC(x, p_dist)
        return np.array([x,p_dist])
        # self.PlotMultipleSim( y, y, np.array([lab,lab]), x_label, y_label, title_string, file_name)
    # # def SolveMultiModel(self,D_p_arr)
    
    # 
    def ParameterEffect(self,scales,D_p=None,V_p = None,half_life = None,Jin = None):
        output = {};
        for ind,scale in enumerate(scales):
            current_D_p = self.D_p
            current_V_p = self.V_p 
            current_half_life = self.half_life
            current_Jin = self.Jin
            # output[ind] = SingleProteinModelWithoutSpines(self.D_p,self.V_p,self.half_life,self.Jin)
            if D_p:
                current_D_p *= scale
            if V_p:
                current_V_p *= scale
            if half_life:
                current_half_life *= scale
            if Jin:
                current_Jin *= scale
            output[ind] = SingleProteinModelWithoutSpines(current_D_p,current_V_p,current_half_life,current_Jin)
            print(output[ind].__dict__)
        return output
    
    def IntegralBC(self,x,p):
        total_p = self.Jin/self.Lamda_p;# + self.Jsin/self.Lamda_ps
        # x,p= self.solveModel()
        num_total_p = (np.sum(p))*(x[1]-x[0])
        print("total p analytic/total p numeric = ",total_p/num_total_p)
    
if __name__ == '__main__':
    SP_model1 = SingleProteinModelWithoutSpines(0.45,0.1,4.35,0.1);
    sim_id = "001";
    x,p_dist = SP_model1.solveModel()
    SP_model1.updateModelParams(V_p=0.001);
    SP_model1.solveModel()
    title_string = "Steady-state spatial distribution \n parameters: $D_p = {%.2f}, V_p = {%.1e}$, half-life = %.2f, Jin= %.2f" \
        %( SP_model1.D_p, SP_model1.V_p, SP_model1.half_life,SP_model1.Jin);
    lab =  'AMPA-R'
    x_label = r'Dendritic distance (in $\mu$M)';
    y_label= r'Protein number';
    folder= "Figures/OneProtein/NoUptake/";
    file_name = folder+"SingleProtein_SingleSim_{0}".format(sim_id);
    pwa = PlottingWidgetAMPA()
    pwa.CreateFolderRecursive(folder)
    norm_types = ["Max","Sum"];
    for norm_type in norm_types:
        pwa.PlotSingleSimSingleProtein(x, p_dist, lab, x_label, y_label, title_string, file_name)
        pwa.plotNormSingle(x, p_dist, lab, x_label, y_label, title_string, file_name,norm_type=norm_type)
        
        scales = np.array([1/5,1/2,1,2])
        modified_objs = SP_model1.ParameterEffect(scales,V_p = 1)
        labels =  ["Vp X {%.2f}"%(scale) for scale in scales]
        title_string = "Effect of changing velocity of active transport";
        all_op = np.zeros((len(scales),len(x)))
        for ind,key in enumerate(modified_objs.keys()):
            x,all_op[ind] = modified_objs[key].solveModel()
        print(all_op)
        file_name = folder+"SingleProtein_MultiSim_velocity_{0}".format(sim_id)
        pwa.PlotMultipleSim(x, all_op,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMulti(x, all_op,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        
        scales = np.array([1/100,1,100,1000])
        title_string = "Effect of changing diffusion constant";
        modified_objs = SP_model1.ParameterEffect(scales,D_p = 1)
        labels =  ["Dp X {%.2f}"%(scale) for scale in scales]
        all_op = np.zeros((len(scales),len(x)))
        for ind,key in enumerate(modified_objs.keys()):
            x,all_op[ind] = modified_objs[key].solveModel()
        print(all_op)
        file_name = folder+"SingleProtein_MultiSim_Diffusion_{0}".format(sim_id)
        pwa.PlotMultipleSim(x, all_op,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMulti(x, all_op,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)
        
        scales = np.array([1/5,1/2,1,2])
        title_string = "Effect of changing half-life";
        modified_objs = SP_model1.ParameterEffect(scales,half_life = 1)
        labels =  ["T-half X {%.2f}"%(scale) for scale in scales]
        all_op = np.zeros((len(scales),len(x)))
        for ind,key in enumerate(modified_objs.keys()):
            x,all_op[ind] = modified_objs[key].solveModel()
        print(all_op)
        file_name = folder+"SingleProtein_MultiSim_half_life_{0}".format(sim_id)
        pwa.PlotMultipleSim(x, all_op,labels,x_label,y_label,title_string,file_name,save_it=1)
        pwa.plotNormMulti(x, all_op,labels,x_label,y_label,title_string,file_name,norm_type=norm_type)