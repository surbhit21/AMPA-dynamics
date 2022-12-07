#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 11:34:03 2022

@author: surbhitwagle
"""


from datetime import datetime
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import RK45
# np.set_printoptions(precision=10)



# simulation parameters
L = 20.0
dx = 0.5
x_grid = np.arange(0,L+dx,dx)
x_points = x_grid.shape[0]

T = 200.0
dt =0.1
t_grid =  np.arange(0,T+dt,dt)
t_points = t_grid.shape[0]


# save after every x percent job is done
job_percent = 10;
t_job_percent = int(T*job_percent/100)
breakpoint()
# model parameters such as diffusion and drift coefficient 

D_s = 0.021   # in uM^2/s
D_c = 0.2   # in uM^2/s
V_p = 4.06e-04    # in uM/s
half_life_surf = np.Inf # in days
Lamda_ps = np.log(2)/(half_life_surf*24*60*60);
half_life_int = 1.95 # in days
Lamda_pc = np.log(2)/(half_life_int*24*60*60);
alpha = 1.0;
beta = alpha/1.15;
Jcin = 10;
Jsin = Jcin/1.15;
omega = 60;             #same as eta_s_max in steady state model file
eta = 1/(15*60);        #same as eta_s_zero in steady state model file
gamma = 1/(50*60);


# cranc-nicolson parameters K_1, K_2, R

K_1 = (D_c*dt)/(2*dx*dx)
K_2 = (D_s*dt)/(2*dx*dx)
R = (V_p*dt)/(4*dx)


# Initial values of P_c, P_s, P_spine is set to 0

P_s_init = np.zeros(x_grid.shape)
P_c_init = np.zeros(x_grid.shape)
P_spine_init = np.zeros(x_grid.shape)

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
    v = u*dt*param
    v[0] += addn
    return v

# implments a single time step integration to get the value of F(P_s,P_{spine})
# def Spine(p_spine,p_s,t,dt):
#     P_spine_next = eta*p_s*(omega) - gamma* p_spine
#     # P_spine_next = -1.0 * P_s_next
    
#     return P_spine_next

def Spine_rk(p_spine, p_s,dt):
    f_vec = lambda p_spine, p_s: eta*p_s*(omega - p_spine) - gamma* p_spine
    k1 = f_vec(p_spine, p_s)
    k2 = f_vec(p_spine + np.multiply(dt/2., k1), p_s - np.multiply(dt/2., k1))
    k3 = f_vec(p_spine + np.multiply(dt/2., k2), p_s - np.multiply(dt/2., k2))
    k4 = f_vec(p_spine + np.multiply(dt, k3), p_s - np.multiply(dt, k3))
    
    return np.multiply(dt/6., k1 + np.multiply(2., k2) + np.multiply(2., k3) + k4)

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
    p_c = np.asarray(pc)
    with open('{0}/PC_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_c)
    f.close()
    p_s = np.asarray(ps)
    with open('{0}/PS_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_s)
    f.close()
    p_sp = np.asarray(pspine)
    with open('{0}/PSPINE_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_sp)
    f.close()
ti=1

for ti in range(1,t_points):
    if (ti*dt/T)*100 % t_job_percent == 0:
        # breakpoint()
        print(ti)
        saveoutput(op_dir,date_time,P_C,P_S,P_SPINE,int((ti*dt/T)*100))
    p_spine_new = Spine_rk(P_spine_now,P_s_now,dt)
    P_c_new = np.linalg.solve(A, B.dot(P_c_now) + BC(P_s_now,dt,beta,2.0*K_1*dx*Jcin/D_c))
    P_s_new =  np.linalg.solve(C, D.dot(P_s_now) + BC(P_c_now,dt,alpha,2.0*K_2*dx*Jsin/D_s)-p_spine_new)
    
    P_c_now = P_c_new
    P_s_now = P_s_new
    P_spine_now = p_spine_new
    
    P_C.append(P_c_now)
    P_S.append(P_s_now)
    P_SPINE.append(P_spine_now)


saveoutput(op_dir,date_time,P_C,P_S,P_SPINE,int((ti*dt/T)*100))
P_C = np.asarray(P_C)
P_S = np.asarray(P_S)
P_SPINE = np.asarray(P_SPINE)

# breakpoint()

# breakpoint()
# plt.ylim((0., 2.1))
plt.xlabel('x'); 
plt.ylabel('concentration')
plt.plot(x_grid, P_c_now,label=r"P_c")
plt.plot(x_grid, P_s_now,label=r"P_s")
plt.plot(x_grid, P_spine_now,label=r"P_{spine}")
plt.legend()
plt.show()

# breakpoint()

fig, (ax0,ax1,ax2) = plt.subplots(1,3)
plt.xlabel('x'); 
plt.ylabel('t')
heatmap1 = ax0.pcolor(x_grid, t_grid, P_C, vmin=P_C.min(), vmax=P_C.max())
heatmap2 = ax1.pcolor(x_grid, t_grid, P_S, vmin=P_S.min(), vmax=P_S.max())
heatmap3 = ax2.pcolor(x_grid, t_grid, P_SPINE, vmin=P_SPINE.min(), vmax=P_SPINE.max())
colorbar1 = plt.colorbar(heatmap1)
colorbar2 = plt.colorbar(heatmap2)
colorbar3 = plt.colorbar(heatmap3)
colorbar1.set_label('Number in cytoplasm')
colorbar2.set_label('Number in surface')
colorbar3.set_label('Number in spine')
plt.show()


