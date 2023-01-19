#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 11:34:03 2022

@author: surbhitwagle
"""

from AMPA_model import RunSim5, RunModelWithFile
from datetime import datetime
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import RK45
np.set_printoptions(precision=5)
# np.set_printoptions(suppress=True)




T = 2.0
dt = 0.01
t_grid =  np.arange(0,T+dt,dt)
t_points = t_grid.shape[0]


# save after every x percent job is done
job_percent = 10;
t_job_percent = int(T*job_percent/100)




# Initial values of P_c, P_s, P_spine is set to 0

P_s_init,P_c_init,P_spine_init,SP_model1 = RunModelWithFile("./ModelParams.json")

# P_s_init = np.zeros(x_grid.shape)
# P_c_init = np.zeros(x_grid.shape)
# P_spine_init = np.zeros(x_grid.shape)

# model parameters such as diffusion and drift coefficient 

D_s = SP_model1.D_s   # in uM^2/s
D_c = SP_model1.D_c  # in uM^2/s
V_p = SP_model1.V_p    # in uM/s
half_life_surf = SP_model1.half_life_surf # in days
Lamda_ps = np.log(2)/(half_life_surf*24*60*60);
half_life_int =  SP_model1.half_life_int # in days
Lamda_pc = np.log(2)/(half_life_int*24*60*60);
alpha =  SP_model1.alpha
beta =  SP_model1.beta;
Jcin =  SP_model1.Jcin;
Jsin =  SP_model1.Jsin;
omega =  SP_model1.omega;             #same as eta_s_max in steady state model file
eta =  SP_model1.eta;        #same as eta_s_zero in steady state model file
gamma =  SP_model1.gamma


# simulation parameters
L = 500.0
dx = SP_model1.dx
x_grid = np.arange(0,L,dx)
x_points = x_grid.shape[0]
# breakpoint()
 # Jcin=0.021
 #    alpha=1
 #    beta = 0.8673779188794962
 #    eta_s0=1/(15*60)
 #    gamma=1/(15*60)
# cranc-nicolson parameters K_1, K_2, R

K_1 = (D_c*dt)/(2*dx*dx)
K_2 = (D_s*dt)/(2*dx*dx)
R = (V_p*dt)/(4*dx)

# matrix parameters
a_1 = (K_1-R)
a_2 = 2.0*K_1 + dt*(Lamda_pc+beta)/2.0
c_1 = K_2
c_2 = 2.0*K_2 + dt*(Lamda_ps+alpha)/2.0

a_2_bc  = (dt*(Lamda_pc+beta)/2.0 + R*(1.0 + V_p*dx/D_c) + K_1 *(1-V_p*dx/D_c))
c_2_bc = K_2 + dt*(Lamda_ps+alpha)/2.0

# creating matrices

A = np.diagflat([-a_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.+a_2_bc]+[1.+a_2 for i in range(x_points-2)]+[1.+a_2_bc]) +\
      np.diagflat([-a_1 for i in range(x_points-1)], 1)
# print(A)
B = np.diagflat([a_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.-a_2_bc]+[1.-a_2 for i in range(x_points-2)]+[1.-a_2_bc]) +\
      np.diagflat([a_1 for i in range(x_points-1)], 1)
      
C =  np.diagflat([-c_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.+c_2_bc]+[1.+c_2 for i in range(x_points-2)]+[1.+c_2_bc]) +\
      np.diagflat([-c_1 for i in range(x_points-1)], 1)

D =  np.diagflat([c_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.-c_2_bc]+[1.-c_2 for i in range(x_points-2)]+[1.-c_2_bc]) +\
      np.diagflat([c_1 for i in range(x_points-1)], 1)
        
# breakpoint()

def BC(u,dt,param,addn):
    f_u = lambda u,param: param*u
    k1 = f_u(u,param)
    k2 = f_u(u+k1*dt/2.0,param)
    k3 = f_u(u+k2*dt/2.0,param)
    k4 = f_u(u+k3*dt,param)
    v = (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
    v[0] += addn;
    return v

# implments a single time step integration to get the value of F(P_s,P_{spine})
# def Spine(p_spine,p_s,t,dt):
#     P_spine_next = eta*p_s*(omega) - gamma* p_spine
#     # P_spine_next = -1.0 * P_s_next
    
#     return P_spine_next

def Spine_rk(p_s, p_spine,dt):
    f_vec = lambda p_s,p_spine: eta*p_s*(omega - p_spine) - gamma* p_spine
    # print(f_vec)
    k1 = f_vec(p_s, p_spine)
    k2 = f_vec(p_s + k1*dt/2.0, p_spine + k1*dt/2.0)
    k3 = f_vec(p_s + k2*dt/2.0, p_spine + k2*dt/2.0)
    k4 = f_vec(p_s + dt*k3, p_spine + k3*dt)
    
    return (dt/6)* (k1 + 2.0*k2 + 2.0*k3 + k4)

P_C = []
P_S = []
P_SPINE = []

P_c_now = P_c_init
P_s_now = P_s_init
P_spine_now = P_spine_init

P_C.append(P_c_now)
P_S.append(P_s_now)
P_SPINE.append(P_spine_now)

#  output folder and file setup
now = datetime.now()
date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
op_dir = os.getcwd()+"/"+date_time
os.makedirs(op_dir, exist_ok=True)
print("date and time:",date_time)	

def saveoutput(op_dir,date_time,pc,ps,pspine,percent):
    p_c = np.asarray(pc,dtype=np.int32)
    with open('{0}/PC_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_c)
    f.close()
    p_s = np.asarray(ps,dtype=np.int32)
    with open('{0}/PS_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_s)
    f.close()
    p_sp = np.asarray(pspine,dtype=np.int32)
    with open('{0}/PSPINE_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_sp)
    f.close()
ti=1

for ti in range(1,t_points):
    # if (ti*dt/T)*100 % t_job_percent == 0:
    #     # breakpoint()
    #     print(ti)
    #     saveoutput(op_dir,date_time,P_C,P_S,P_SPINE,int((ti*dt/T)*100))
    delta_pspine = Spine_rk(P_s_now,P_spine_now,dt)
    # breakpoint()
    # print(delta_pspine[:10])
    
    P_spine_now = P_spine_now + delta_pspine
    
    P_c_new = np.linalg.solve(A, B.dot(P_c_now) + BC(P_s_now,dt,alpha,2.0*K_1*dx*Jcin/D_c))
    P_s_new =  np.linalg.solve(C, (D.dot(P_s_now) + BC(P_c_now,dt,beta,2.0*K_2*dx*Jsin/D_s) - delta_pspine))
    # print(((P_c_new-P_c_now)/(P_c_now))[:10])
    
    # print(((P_s_new-P_s_now)/(P_s_now))[:10])
    # print("*"*10)
    # print(sum(p_spine_new/P_spine_init),sum(P_s_new/P_s_init),sum(P_c_new/P_c_init))
    P_c_now = P_c_new
    P_s_now = P_s_new
    # P_spine_now = P_spine_now
    
    P_C.append(P_c_now)
    P_S.append(P_s_now)
    P_SPINE.append(P_spine_now)



saveoutput(op_dir,date_time,P_C,P_S,P_SPINE,100)
# breakpoint()



