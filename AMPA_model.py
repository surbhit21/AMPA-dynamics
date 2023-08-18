#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 16:22:34 2022

@author: surbhitwagle
using code styling from PEP8 - https://peps.python.org/pep-0008/

the main file used to generate mRNA plots for the AMPA paper
"""

import json
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import math
import numpy as np
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
from scipy.integrate import solve_bvp
import scipy.optimize as opt       # import root-finding algorithm
# import sympy as sp                 # Python toolbox for symbolic maths



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

def PlotSingleSimSingleProtein(x, y1, y2, lab1, lab2,
                               xlab, ylab, title_string,
                               file_name, width=10, height=8, fsize=16, save_it = 1):
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


def PlotSingleSimTwoProtein(x, ps, pc, p_spine, lab_ps, lab_pc, lab_p_spine,
                            xlab, ylab, title_string,
                            file_name, width=10, height=8, fsize=16, save_it = 1):
    fig, axes = plt.subplots(1, 2,figsize=(width, height))
    ax = axes.ravel()
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
    """
        Class for the steady state solution of  model equations:
            \frac{\partial P_c}{dt} = D_c \frac{\partial^2 P_c}{\partial x^2} - V_p\frac{\partial P_c}{\partial x} 
            - \lambda_c P_c - \beta P_c + \alpha P_s
            
            \frac{\partial P_s}{dt} &= D_s \frac{\partial^2 P_s}{\partial x^2} - \lambda_s P_s -\alpha P_s 
            + \beta P_c -\eta P_s (\omega- P_{spine}) + \gamma P_{spine}  
             
            \frac{\partial P_{spine}}{dt} &= \eta P_s (\omega- P_{spine}) - \gamma P_{spine} 
    """

    def __init__(self, D_s, D_c, V_p, half_life_surf, half_life_int, alpha,
                 beta, Jsin, Jcin, omega, eta,gamma,dx):
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
        self.omega = omega;
        self.eta = eta;
        self.gamma = gamma;
        self.dx = dx
        self.x_grid = np.arange(0,L,dx)
        
    
    def fun(self, x, y):
        """
        function that deines the darivatives
        """
        ps,dps,pc,dpc = y
        return [dps,\
                ((self.alpha+self.Lamda_ps)/self.D_s)*ps  - (self.beta/self.D_s)*pc, dpc,\
                ((self.beta+self.Lamda_pc)/self.D_c)*pc + (self.V_p/self.D_c)*dpc - (self.alpha/self.D_c)*ps\
               ]
                             
        
    def bc(self,ya,yb):
        """
            function for boundary condition values
        """
        return np.array([ya[1] + self.Jsin/self.D_s, yb[1],  ya[3] - self.V_p*ya[2]/self.D_c + self.Jcin/self.D_c,yb[3]- self.V_p*yb[2]/self.D_c])

    def solveNumerical(self):
        y = np.zeros((4,self.x_grid.size))
        soln = solve_bvp(self.fun, self.bc, self.x_grid, y,max_nodes=1e+9, verbose=1,tol=1e-3, bc_tol=1e-8)
        # breakpoint()
        ps_dist = soln.y[0]
        pc_dist = soln.y[2]
        self.x_grid = soln.x
        # breakpoint()
        print(self.bc(soln.y[:,0],soln.y[:,-1]))
        return ps_dist, pc_dist
        
    def GetODERoots(self):
        A =  [[0, 0, 1, 0],
              [0, 0, 0, 1],
              [(self.alpha+self.Lamda_ps)/self.D_s, -self.beta/self.D_s, 0, 0],
              [-self.alpha/self.D_c, (self.beta+self.Lamda_pc)/self.D_c, 0, self.V_p/self.D_c]];
        # breakpoint()
        ss_eigen_values,ss_eig_vec = np.linalg.eig(A);
        return ss_eigen_values,ss_eig_vec
   
    def Descriminant(self):
        a = 1;
        b = -self.V_p/self.D_c;
        c = -((self.beta+self.Lamda_pc)/self.D_c + (self.alpha+self.Lamda_ps)/self.D_s);
        d = self.V_p*(self.alpha+self.Lamda_ps)/(self.D_c*self.D_s);
        e = (self.Lamda_ps*self.beta + self.alpha*self.Lamda_pc + self.Lamda_pc*self.Lamda_ps)/(self.D_s*self.D_c);
        
        des = 256*(a**3)*(e**3) - 192*(a**2)*(b*d)*(e**2) - 128*(a*c*e)**2 + 144*(a*d)**2*(c*e) - 27*(a**2)*(d**4) 
        + 144*a*c*((b*e)**2) - 6*a*e*((b*d)**2) - 80*a*b*d*e*(c**2) + 18*a*b*c*(d**3) + 16*a*e*(c**4) - 4*a*(c**3)*(d**2) 
        - 27*(b**4)*(e**2) + 18*c*d*e*(b**3) - 4*(b*d)**3 - 4*e*(b**2)*(c**3) + (b*c*d)**2
        
        print("Descriminant = ",des)
        return des
    
    def GetODEPrefectors(self,roots):
        r1,r2,r3,r4 = roots
        A = np.array([[self.D_s*r1**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r1**2 - self.V_p*r1 - self.Lamda_pc, 0, 0, 0],
             [0,self.D_s*r2**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r2**2 - self.V_p*r2 - self.Lamda_pc, 0, 0, ],
             [0, 0, self.D_s*r3**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r3**2 - self.V_p*r3 - self.Lamda_pc, 0 ],
             [0, 0, 0, self.D_s*r4**2 - self.Lamda_ps, 0, 0, 0, self.D_c*r4**2 - self.V_p*r4 - self.Lamda_pc],
             [r1*np.exp(r1*L), r2*np.exp(r2*L), r3*np.exp(r3*L), r4*np.exp(r4*L), 0, 0, 0, 0],
             [r1, r2, r3, r4, 0, 0, 0, 0],
             [0, 0, 0, 0,(self.D_c*r1 - self.V_p)*np.exp(r1*L), (self.D_c*r2 - self.V_p)*np.exp(r2*L), (self.D_c*r3 - self.V_p)*np.exp(r3*L), (self.D_c*r4 - self.V_p)*np.exp(r4*L) ],
             [0, 0, 0, 0,self.D_c*r1 - self.V_p, self.D_c*r2 - self.V_p, self.D_c*r3 - self.V_p, self.D_c*r4 - self.V_p ]]);
        B = np.array([0, 0, 0, 0, 0, -self.Jsin, 0, -self.Jcin]);
        # breakpoint()
        prefactors = np.linalg.solve(A,B)
        return prefactors
    
    def SolveAnalytical(self):
        desc = self.Descriminant()
        roots,vecs = self.GetODERoots()
        pfs = self.GetODEPrefectors(roots)
        ps_dist = np.zeros((np.shape(self.x_grid)))
        pc_dist = np.zeros((np.shape(self.x_grid)))
        for i in range(0,4):
            print(i)
            ps_dist += pfs[i]*np.exp(roots[i]*self.x_grid)
            pc_dist += pfs[4+i]*np.exp(roots[i]*self.x_grid)
        return ps_dist,pc_dist
    
    def IntegralBC(self,ps_dist,pc_dist):
        tot_ps_num = np.sum(ps_dist)*self.dx
        tot_pc_num = np.sum(pc_dist)*self.dx
        tot_ps_ana = (self.beta*self.Jcin)/(self.Lamda_pc*self.alpha)
        tot_pc_ana = self.Jcin/self.Lamda_pc
        return tot_ps_ana/tot_ps_num,tot_pc_ana/tot_pc_num
    
def SaveFigures(filename,ext_list = [".png",".svg",".pdf"]):
    for ext in ext_list:
        plt.savefig(filename+ext,dpi=300)
        

def RunSim5(delta_x,v_p,D_c,D_s):
    Jcin= 0.05
    alpha= 1.5e-4
    beta = alpha*2
    eta_s0= 1e-3
    gamma= 1/(43) 
    SP_model1 = DendriteWithStochasticSpinesConstantV(D_s,D_c,v_p,float('inf'),4.5,alpha,beta,0,Jcin,60,eta_s0,gamma,delta_x);
    sim_id = "002";
    ps_dist,pc_dist = SP_model1.solveNumerical()
    # x=np.arange(0,L,delta_x)
    SP_model1.x_grid = SP_model1.x_grid.tolist()
    # jsonstr1 = json.dumps(SP_model1.__dict__)
    param_dict = {}
    param_dict["D_c"] = SP_model1.D_c
    param_dict["D_s"] = SP_model1.D_s
    param_dict["V_p"] = SP_model1.V_p
    param_dict["half_life_surf"] = SP_model1.half_life_surf
    param_dict["half_life_int"] = SP_model1.half_life_int
    param_dict["alpha"] = SP_model1.alpha
    param_dict["beta"] = SP_model1.beta
    param_dict["gamma"] = SP_model1.gamma
    param_dict["eta"] = SP_model1.eta
    param_dict["omega"] = SP_model1.omega
    param_dict["dx"] = SP_model1.dx
    param_dict["Jsin"] = SP_model1.Jsin
    param_dict["Jcin"] = SP_model1.Jcin
    # breakpoint()
    with open ("./ModelParams.json","w") as fp:
        json.dump(param_dict,fp)
    fp.close()
    ps_spine = SP_model1.omega*(1/(1+ (SP_model1.gamma/(SP_model1.eta*ps_dist))))
    # fig,ax = plt.subplots(figsize=(10, 8))
    # fsize=16
    # ax.plot(SP_model1.x_grid,ps_dist,label=r"$p_s$",color = color_surf,linewidth=3.0)
    # ax.plot(SP_model1.x_grid,pc_dist,label=r"$P_c$",color = color_cyto,linewidth=3.0)
    # ax.plot(SP_model1.x_grid,ps_spine,label=r"$P_{spine}$",color = color_spine,linewidth=3.0)
    # fig.tight_layout()
    # plt.legend(prop={'size': fsize})
    # SaveFigures("./ModelDist")
    # plt.show()
    # breakpoint()
    return ps_dist,pc_dist,ps_spine

RunSim5(0.24,4.06e-04,0.2,0.2)


def RunModelWithFile(param_file):
     with open (param_file,"r") as fp:
        params = json.load(fp)
     fp.close()
     SP_model1 = DendriteWithStochasticSpinesConstantV(**params);
     SP_model1.x_grid = np.asarray(SP_model1.x_grid)
     ps_dist,pc_dist = SP_model1.solveNumerical()
     ps_spine = SP_model1.omega*(1/(1+ (SP_model1.gamma/(SP_model1.eta*ps_dist))))
     return ps_dist, pc_dist, ps_spine, SP_model1