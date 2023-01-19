#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 11:34:03 2022

@author: surbhitwagle
"""

from AMPA_model import RunSim5
from datetime import datetime
import os
import matplotlib.pyplot as plt
import numpy as np
# np.set_printoptions(precision=10)



simulation parameters
L = 500.0
dx = 0.1
x_grid = np.arange(0,L,dx)
x_points = x_grid.shape[0]

T = 50
dt =0.1
t_grid =  np.arange(0,T+dt,dt)
t_points = t_grid.shape[0]

# save after every x percent job is done
job_percent = 10;
t_job_percent = int(T*job_percent/100)

# model parameters such as diffusion and drift coefficient 

D_s = 0.021   # in uM^2/s
D_c = 0.2   # in uM^2/s
V_p = 4.06e-04    # in uM/s
half_life_surf = np.Inf # in days
Lamda_ps = np.log(2)/(half_life_surf*24*60*60);
half_life_int = 1.95 # in days
Lamda_pc = np.log(2)/(half_life_int*24*60*60);
alpha = 1.0;
beta = 0.8673779188794962;
Jcin = 0.021;
Jsin = 0.01637936;
omega = 0;             #same as eta_s_max in steady state model file
eta = 1/(15*60);        #same as eta_s_zero in steady state model file
gamma = 1/(15*60);




# cranc-nicolson parameters K_1, K_2, R

K_1 = (D_c*dt)/(2*dx*dx)
K_2 = (D_s*dt)/(2*dx*dx)
R = (V_p*dt)/(4*dx)


# Initial values of P_s and P_c is set to 0
P_s_init,P_c_init,P_spine_init = RunSim5(dx,V_p,D_c,D_s)
# P_s_init = np.zeros(x_grid.shape)
# P_c_init = np.zeros(x_grid.shape)

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

def F(u,dt,param,addn):
    v = u*dt*param
    v[0] += addn
    return v


P_C = []
P_S = []

P_c_now = P_c_init
P_s_now = P_s_init

P_C.append(P_c_now)
P_S.append(P_s_now)

#  output folder and file setup
now = datetime.now()
date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
op_dir = os.getcwd()+"/"+date_time
# breakpoint()
os.makedirs(op_dir, exist_ok=True)
print("date and time:",date_time)	

def saveoutput(op_dir,date_time,pc,ps,percent):
    p_c = np.asarray(pc)
    with open('{0}/PC_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_c)
    f.close()
    p_s = np.asarray(ps)
    with open('{0}/PS_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_s)
    f.close()
   
ti=1


for ti in range(1,t_points):
    # if (ti*dt/T)*100 % t_job_percent == 0:
    #     # breakpoint()
    #     # print(ti)
    #     saveoutput(op_dir,date_time,P_C,P_S,int((ti*dt/T)*100))
    P_c_new = np.linalg.solve(A, B.dot(P_c_now) + F(P_s_now,dt,beta,2.0*K_1*dx*Jcin/D_c))
    P_s_new =  np.linalg.solve(C, D.dot(P_s_now) + F(P_c_now,dt,alpha,2.0*K_2*dx*Jsin/D_s))
    
    P_c_now = P_c_new
    P_s_now = P_s_new

    
    P_C.append(P_c_now)
    P_S.append(P_s_now)



saveoutput(op_dir,date_time,P_C,P_S,100)

P_C = np.asarray(P_C)
P_S = np.asarray(P_S)


# plt.ylim((0., 2.1))
plt.xlabel('x'); 
plt.ylabel('concentration')
plt.plot(x_grid, P_c_init,label="pc0")
plt.plot(x_grid, P_s_init,label="ps0")
plt.plot(x_grid, P_c_now,label="pct")
plt.plot(x_grid, P_s_now,label="pst")
plt.legend()
plt.show()

# breakpoint()

fig, (ax0,ax1) = plt.subplots(1,2)
plt.xlabel('x'); 
plt.ylabel('t')
heatmap = ax1.pcolor(x_grid, t_grid, P_C, vmin=P_C.min(), vmax=P_C.max())
heatmap = ax0.pcolor(x_grid, t_grid, P_S, vmin=P_S.min(), vmax=P_S.max())
colorbar = plt.colorbar(heatmap)
colorbar.set_label('concentration U')
plt.show()


