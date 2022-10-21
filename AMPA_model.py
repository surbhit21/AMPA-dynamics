#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 16:22:34 2022

@author: surbhitwagle
"""

from scipy.integrate import solve_bvp
import numpy as np
import scipy.optimize as opt       # import root-finding algorithm
import sympy as sp                 # Python toolbox for symbolic maths
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D # Toolbox for rendring 3D figures
# from mpl_toolkits import mplot3d   # Toolbox for rendring 3D figures
# @title Figure Settings
# import ipywidgets as widgets  # interactive display
# from ipywidgets import interact
# from pathlib import Path
import math
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
# %config InlineBackend.figure_format = 'retina'
# use NMA plot style
# plt.style.use("https://raw.githubusercontent.com/NeuromatchAcademy/course-content/master/nma.mplstyle")
# my_layout = widgets.Layout()

fig_w, fig_h = 12, 4.5
my_fontsize = 16
my_params = {'axes.labelsize': my_fontsize,
          'axes.titlesize': my_fontsize,
          'figure.figsize': [fig_w, fig_h],
          'font.size': my_fontsize,
          'legend.fontsize': my_fontsize-4,
          'lines.markersize': 8.,
          'lines.linewidth': 2.,
          'xtick.labelsize': my_fontsize-2,
          'ytick.labelsize': my_fontsize-2}

plt.rcParams.update(my_params)

color_surf = '#005f73'
color_cyto = '#9b2226'
color_spine = '#CA6702'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
color_list = [color_surf, color_cyto, color_spine, CB91_Amber, CB91_Purple, CB91_Violet]

def PlotSingleSimSingleProtein(x,y1,y2,lab1,lab2,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
    print(width,height)
    fig, axes = plt.subplots(1, 2,figsize=(width, height))
    ax = axes.ravel()
    fig.suptitle(title_string,fontsize=fsize)
#     plt.rc('font', **{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=False)
    ax[0].plot(x,y1,label=lab1)
    ax[0].set_xlabel(xlab,fontsize=fsize)
    ax[0].set_ylabel(ylab,fontsize=fsize)
#     plt.title(title_string,fontsize=fsize)
    ax[1].plot(x,y2,label=lab2)
    ax[1].set_xlabel(xlab,fontsize=fsize)
    ax[1].set_ylabel(ylab,fontsize=fsize)
    
    plt.legend(prop={'size': fsize})
    # plt.show()
    folder = "."
    if save_it == 1:
        plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=150)
        plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=150)
        print("saved figures to: {%s/%s}" %(folder, file_name))
    else:
        print("Plots not saved")
    plt.show()


def PlotSingleSimTwoProtein(x,ps,pc,p_spine,lab_ps,lab_pc,lab_p_spine,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
    fig, axes = plt.subplots(1, 2,figsize=(width, height))
    ax = axes.ravel()
#     plt.rc('font', **{'family':'serif','serif':['Palatino']})
#     plt.rc('text', usetex=True)
    fig.suptitle(title_string,fontsize=fsize)
    ax[0].plot(x,ps,label=lab_ps,color=color_list[0])
    ax[0].plot(x,pc,label=lab_pc,color=color_list[1])
    ax[0].set_xlabel(xlab,fontsize=fsize)
    ax[0].set_ylabel(ylab,fontsize=fsize)
#     plt.title(title_string,fontsize=fsize)
    ax[0].legend(prop={'size': fsize})
#     ax[0].ylim(0,3)
#     plt.ylim(0,1.1)
    ax[1].plot(x,p_spine,label=lab_p_spine,color=color_list[2])
    ax[1].set_xlabel(xlab,fontsize=fsize)
    ax[1].set_ylabel('Filling fraction',fontsize=fsize)
    
#     plt.title(title_string,fontsize=fsize)
    ax[1].legend(prop={'size': fsize})
#     plt.ylim(0,2)
    # plt.show()
    folder = "."
    if save_it == 1:
        plt.savefig("%s/%s.%s"%(folder,file_name,"png"),dpi=150)
        plt.savefig("%s/%s.%s"%(folder,file_name,"eps"),dpi=150)
        print("saved figures to: {%s/%s}" %(folder, file_name))
    else:
        print("Plots not saved")
    plt.show()
    

def CreateFolderRecursive(folder):
        Path(folder).mkdir(parents=True, exist_ok=True)

L= 500
class DendriteWithStochasticSpinesConstantV():
    def __init__(self,D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin,eta_s_max,eta_s_zero,gamma):
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
        self.eta_s_max = eta_s_max;
        self.eta_s_zero = eta_s_zero;
        self.gamma = gamma;
        
    def updateModelParams(self,D_s = None,D_c = None,V_p = None,half_life_surf = None,half_life_int = None\
                          ,alpha = None,beta = None,Jsin = None,Jcin = None,eta_s_max=None,eta_s_zero = None,gamma=None):
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
        if eta_s_max:
            self.eta_s_max = eta_s_max;
#         if eta_c_max:
#             self.eta_c_max = eta_c_max;
        if eta_s_zero:
            self.eta_s_zero = eta_s_zero;
        if gamma:
            self.gamma = gamma;
    
    def fun(self,x,y):
        ps,dps,pc,dpc = y
        # print(dps[0],dps[-1],dpc[0],dpc[-1])
        return [dps,\
                ((self.alpha+self.Lamda_ps)/self.D_s)*ps  - (self.beta/self.D_s)*pc, dpc,\
                ((self.beta+self.Lamda_pc)/self.D_c)*pc + (self.V_p/self.D_c)*dpc - (self.alpha/self.D_c)*ps\
               ]
                             
        
    def bc(self,ya,yb):
        return np.array([self.D_s*ya[1] + self.Jsin,  self.D_c*ya[3] - self.V_p*ya[2] + self.Jcin, self.D_s*yb[1],self.D_c*yb[3]- self.V_p*yb[2]])

    def solveModel(self,delta_x):
        # delta_x = 0.0114; #step size for the numerical simulation
        # sim_id="001";
        # print(len(x_line))
        # solving model
        # D_p * p'' - (V_p(x)*p)' - Lamda_p*p = 0
        x=np.arange(0,L,delta_x)
#         print(x)
        # params=np.array([L]);
        y = np.zeros((4,x.size))
        soln = solve_bvp(self.fun, self.bc, x, y,max_nodes=1e+9)
#         print(len(soln))
        ps_dist = soln.sol(x)[0]
        pc_dist = soln.sol(x)[2]
#         ps_spine = soln.sol(x)[4]
#         self.IntegralBC(delta_x,ps_dist,pc_dist)
        return x, ps_dist, pc_dist
    
#     def IntegralBC(self,delta_x,ps_dist,pc_dist):
#         protein_loss = -self.Jcin + (self.Lamda_pc*np.sum(pc_dist) + np.sum(self.eta_s_max*np.tanh(self.eta_s_zero*ps_dist)))*delta_x
#         print("protein lost/found  in the void = ",protein_loss)
        
    def GetODERoots(self):
        A =  [[0, 0, 1, 0],
              [0, 0, 0, 1],
              [(self.alpha+self.Lamda_ps)/self.D_s, -self.beta/self.D_s, 0, 0],
              [-self.alpha/self.D_c, (self.beta+self.Lamda_pc)/self.D_c, 0, self.V_p/self.D_c]];
        breakpoint()
        ss_eigen_values,ss_eig_vec = np.linalg.eig(A);
#         print(A)
#         print(ss_eigen_values)
        return ss_eigen_values,ss_eig_vec
    def Descriminant(self):
        a = 1;
        b = -self.V_p/self.D_c;
        c = -((self.beta+self.Lamda_pc)/self.D_c + (self.alpha+self.Lamda_ps)/self.D_s);
        d = self.V_p*(self.alpha+self.Lamda_ps)/(self.D_c*self.D_s);
        e = (self.Lamda_ps*self.beta + self.alpha*self.Lamda_pc + self.Lamda_pc*self.Lamda_ps)/(self.D_s*self.D_c);
        
        des = 256*(a**3)*(e**3) - 192*(a**2)*(b*d)*(e**2) - 128*(a*c*e)**2 + 144*(a*d)**2*(c*e) - 27*(a**2)*(d**4) \
        + 144*a*c*(b*e)**2 - 6*a*e*(b*d)**2 - 80*a*b*d*e*(c**2) + 18*a*b*c*(d**3) + 16*a*e*(c**4) - 4*a*(c**3)*(d**2) \
        -27*(b**4)*(e**2) + 18*c*d*e*(b**3) - 4*(b*d)**3 - 4*e*(b**2)*(c**3) + (b*c*d)**2
        
        print(des)
        return des
    def GetODEPrefectors(self,roots):
        # roots = roots
        r1,r2,r3,r4 = roots
        A = np.array([[self.D_s*r1**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r1**2 - self.V_p*r1 - self.Lamda_pc, 0, 0, 0],
             [0,self.D_s*r2**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r2**2 - self.V_p*r2 - self.Lamda_pc, 0, 0, ],
             [0, 0, self.D_s*r3**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r3**2 - self.V_p*r3 - self.Lamda_pc, 0 ],
             [0, 0, 0, self.D_s*r4**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r4**2 - self.V_p*r4 - self.Lamda_pc],
             [r1*np.exp(r1*L), r2*np.exp(r2*L), r3*np.exp(r3*L), r4*np.exp(r4*L), 0, 0, 0, 0],
             [r1, r2, r3, r4, 0, 0, 0, 0],
             [0, 0, 0, 0,(self.D_c*r1 - self.V_p)*np.exp(r1*L), (self.D_c*r2 - self.V_p)*np.exp(r2*L), (self.D_c*r3 - self.V_p)*np.exp(r3*L), (self.D_c*r4 - self.V_p)*np.exp(r4*L) ],
             [0, 0, 0, 0,self.D_c*r1 - self.V_p, self.D_c*r2 - self.V_p, self.D_c*r3 - self.V_p, self.D_c*r4 - self.V_p ]]);
        B = np.array([0, 0, 0, 0, 0, 0, 0, -self.Jcin]);
        breakpoint()
        prefactors = np.linalg.solve(A,B)
#         print(prefactors)
        return prefactors
    def SolveAnalytical(self,delta_x):
        desc = self.Descriminant()
        roots,vecs = self.GetODERoots()
        # breakpoint()
        pfs = self.GetODEPrefectors(roots)
        # delta_x = 0.0114;
        x=np.arange(0,L,delta_x)
        ps_dist = np.zeros((np.shape(x)))
        pc_dist = np.zeros((np.shape(x)))
        for i in range(0,4):
            print(i)
            ps_dist += pfs[i]*np.exp(roots[i]*x)
            pc_dist += pfs[4+i]*np.exp(roots[i]*x)
#         print("*****\n",ps_dist,pc_dist)
        # self.IntegralBC(ps_dist,pc_dist,delta_x)
        return x,ps_dist,pc_dist
    
    def IntegralBC(self,ps_dist,pc_dist,delta_x):
        tot_ps_num = np.sum(ps_dist)*delta_x
        tot_pc_num = np.sum(pc_dist)*delta_x
        tot_ps_ana = (self.beta*self.Jcin)/(self.Lamda_pc*self.alpha)
        tot_pc_ana = self.Jcin/self.Lamda_pc
        # print("total_p = ",tot_ps_ana+tot_pc_ana,tot_ps_num+tot_pc_num)
        # print("ratio cytoplasm = ",tot_pc_ana/tot_pc_num)
        # print("ratio surf = ",tot_ps_ana/tot_ps_num)
        return tot_ps_ana,tot_pc_ana
def SaveFigures(filename,ext_list = [".png",".svg",".pdf"]):
    for ext in ext_list:
        plt.savefig(filename+ext,dpi=300)
def RunSim5(delta_x,v_p,D_c,D_s):
    # print("=========================================================================================")
    # print(" v_p =%.5f \t|| Jin = %0.2f \t|| alpha =%0.2f \t|| beta =%0.2f \t|| eta_s0 = %.2e \n"%(v_p,Jcin,alpha,beta,eta_s0),end="")
    # print("=========================================================================================")
    # v_p=0
    Jcin=0.021
    alpha=1
    beta = 0.8673779188794962
    eta_s0=1/(15*60)
    gamma=1/(15*60)  #D_s,D_c,V_p,half_life_surf,half_life_int,alpha,beta,Jsin,Jcin,eta_s_max,eta_s_zero,gamma
    SP_model1 = DendriteWithStochasticSpinesConstantV(D_s,D_c,v_p,float('inf'),1.95,alpha,beta,0.01637936,Jcin,60,eta_s0,gamma);
    sim_id = "002";
    x,ps_dist,pc_dist = SP_model1.solveModel(delta_x)
    # breakpoint()
    ps_spine = SP_model1.eta_s_max*(1/(1+ (SP_model1.gamma/(SP_model1.eta_s_zero*ps_dist))))
    fig,ax = plt.subplots(figsize=(10, 8))
    fsize=16
    ax.plot(x,ps_dist/ps_dist[0],label=r"$p_s$",color = color_surf,linewidth=3.0)
    ax.plot(x,pc_dist/pc_dist[0],label=r"$P_c$",color = color_cyto,linewidth=3.0)
    ax.plot(x,ps_spine/SP_model1.eta_s_max,label=r"$P_{spine}$",color = color_spine,linewidth=3.0)
    # ax.set_ylabel("GluA2 count",fontsize=fsize)
    # ax.set_xlabel('Dendritic distance (in uM)',fontsize=fsize)
    fig.tight_layout()
    plt.legend(prop={'size': fsize})
    SaveFigures("./ModelDist")
    plt.show()
    # breakpoint()
    # breakpoint()
# #     title_string = (r"Steady-state spatial distribution"+" \n parameters:\
# #        "+r" $D_s$ = %.2f, half-life-surf = %.2f, $Jsin= %.2f,  \alpha = %.2f, \eta_{smax} = %.1e ,\eta_{s0}$ = %.1e"+" \n"+ \
# #         r"$D_c = {%.2f}, V_p = {%.1e}$, half-life-int = %.2f, $Jcin= %.2f, \beta = %.2f$") \
# #     %( SP_model1.D_s, SP_model1.half_life_surf,SP_model1.Jsin,SP_model1.alpha,SP_model1.eta_s_max,SP_model1.eta_s_zero,\
# #       SP_model1.D_c,SP_model1.V_p, SP_model1.half_life_int,SP_model1.Jcin,SP_model1.beta);
#     lab_ps =  r'$P_s$'
#     lab_pc =  r'$P_c$'
#     x_label = r'Dendritic distance in ($\mu$M)';
#     y_label= r'distribution of $P_s$ and $P_c$';
#     folder= "Figures/TwoProtein/WithCappedUptake/{0}/".format(sim_id);
#     file_name = folder+ "TwoProtein_SingleSim_withCappedSpinesUptake_{0}".format(sim_id);
# #     pwa = PlottingWidgetAMPA()
#     CreateFolderRecursive(folder)
# #     pwa.PlotSingleSimTwoProtein(x, ps_dist,pc_dist, lab_ps,lab_pc, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
# #     title_string = (r"Steady-state spatial distribution in Spines"+" \n parameters:\
# #        "+r" $D_s = %.2f,    \alpha = %.2f, \eta_{smax} = %.1e ,\eta_{s0} = %.1e, \gamma$ = %.1e"+" \n"+ \
# #         r"$D_c = {%.2f}, V_p = {%.1e}, \lambda_c = {%.1e}, Jcin= %.2f, \beta = %.2f$") \
# #     %( SP_model1.D_s,SP_model1.alpha,SP_model1.eta_s_max,SP_model1.eta_s_zero,SP_model1.gamma,\
# #       SP_model1.D_c,SP_model1.V_p, SP_model1.Lamda_pc,SP_model1.Jcin,SP_model1.beta);
#     title_string = "Steady-state spatial distribution"
#     lab_ps_spine =  r'$P_{spine}$'
# #     lab_pc_spine =  'Spine cytoplsmic AMPA-R'
#     x_label = r'Dendritic distance in ($\mu$M)';
#     y_label= r'copy number';
#     # folder= "Figures/OneProtein/WithCappedUptake/";
#     file_name = folder+"Spine_SingleSim_TwoProtein_capped_uptake_dist_constantV_{}_v_{:e}".format(sim_id,v_p);
# #     ps_spine = SP_model1.eta_s_max*np.tanh(SP_model1.eta_s_zero*ps_dist) 
# #     print("ratio of ps/pc = ",np.sum(ps_dist)/np.sum(pc_dist),"\n beta/alpha =",beta/alpha)
#     ps_spine = SP_model1.eta_s_max*(1/(1+SP_model1.gamma/(SP_model1.eta_s_zero*ps_dist)))
#     print(ps_dist[:20],"\n",pc_dist[:20],"\n",ps_spine[:20])
#     # PlotSingleSimTwoProtein(x,ps_dist,pc_dist, ps_spine/SP_model1.eta_s_max, lab_ps,lab_pc, lab_ps_spine, x_label, y_label, title_string, file_name,fsize=14,save_it = 1)
#     print("Soma Normalized")
# #     PlotSingleSimTwoProtein(x,ps_dist/ps_dist[0],pc_dist/ps_dist[0], ps_spine/ps_spine[0], lab_ps,lab_pc, lab_ps_spine, x_label, y_label, title_string, file_name,fsize=14,save_it = 0)

#     y_label= r'Distribution';
# #     PlotSingleSimTwoProtein(x,ps_dist,pc_dist, ps_spine, lab_ps,lab_pc, lab_ps_spine, x_label, y_label, title_string, file_name,fsize=14,save_it = 0)
    # delta_x = 0.0114
#     PlotSingleSimTwoProtein(x[1:int(160/delta_x)],ps_dist[1:int(160/delta_x)],pc_dist[1:int(160/delta_x)], ps_spine[1:int(160/delta_x)], lab_ps,lab_pc, lab_ps_spine, x_label, y_label, title_string, file_name,fsize=14,save_it = 0)
#     print("sum P_spine/Sum P_s = " , ps_spine.sum()/ps_dist.sum())
    return ps_dist,pc_dist
# _ = widgets.interact(RunSim5,v_p = (0,.01,0.001), Jcin = (0.001,0.1,0.005), alpha = (0.001,1,0.1),beta = (0.01,1,0.1),eta_s0 = (0, 0.1, .01),gamma=(1/43,0.1,0.01))
# RunSim5(0.24,0.02,0.052,0.021)
# RunSim5(0.24,4.06e-04,1.0654e-05,0.26144088)
# RunSim5(0.24,4.06e-03,3,0.26144088)
# RunSim5(0.24, 0.17, 0.14 , 0.54)
RunSim5(0.24,0.042,1.5,6.7)
# ps1,pc1 = RunSim5(0.05,0.021,alpha=0.01,beta = 0.02,eta_s0=1/100,gamma=1/43)
# ps_dist,pc_dist = RunSim5(0.02,1e-3,0.021)
# breakpoint()
# RunSim5(1e-4,0.021,alpha=0.01,beta = 0.02,eta_s0=1/100,gamma=1/43)
# RunSim5(0,0.021,alpha=0.01,beta = 0.02,eta_s0=1/100,gamma=1/43)
