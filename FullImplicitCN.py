#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 15:12:43 2022

@author: surbhitwagle
"""


from AMPA_model import RunSimGluA2, RunModelWithFile
from datetime import datetime
# from LoadNPlot import CallWithRunId
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.integrate import RK45
from scipy.sparse.linalg import spsolve
np.set_printoptions(precision=15)
np.set_printoptions(suppress=True)




T = 100.0
dt = 1
t_grid =  np.arange(0,T+dt,dt)
t_points = t_grid.shape[0]


# save after every x percent job is done
job_percent = 10;
t_job_percent = int(T*job_percent/100)




# Initial values of P_c, P_s, P_spine is set to 0

P_s_init,P_c_init,P_spine_init,SP_model1 = RunModelWithFile("./ModelParams.json")

L = 500.0
dx = SP_model1.dx
x_grid = np.arange(0,L,dx)
x_points = x_grid.shape[0]

# P_s_init = np.zeros(x_grid.shape)
# P_c_init = np.zeros(x_grid.shape)
# P_spine_init = np.zeros(x_grid.shape)

# model parameters such as diffusion and drift coefficient 

D_s = SP_model1.D_s   # in uM^2/s
D_c = SP_model1.D_c   # in uM^2/s
V_p = SP_model1.V_p   # in uM/s
half_life_surf = SP_model1.half_life_surf # in days
Lamda_ps = np.log(2)/(half_life_surf*24*60*60);
half_life_int =  SP_model1.half_life_int  # in days
Lamda_pc = np.log(2)/(half_life_int*24*60*60);
alpha =  SP_model1.alpha 
beta =  SP_model1.beta;
Jcin =  SP_model1.Jcin;
Jsin =  SP_model1.Jsin;
omega =  SP_model1.omega;    # same as eta_s_max in steady state model file
eta =  SP_model1.eta;        # same as eta_s_zero in steady state model file
gamma =  SP_model1.gamma


# simulation parameters




# matrix parameters
K_1 = (D_c*dt)/(2*dx*dx)
K_2 = (D_s*dt)/(2*dx*dx)
R = (V_p*dt)/(4*dx)


a_1 = (K_1-R)
a_2 = 2.0*K_1 + dt*(Lamda_pc+beta)/2.0
c_1 = K_2
c_2 = 2.0*K_2 + dt*(Lamda_ps+alpha)/2.0

# """
# BC forward and backward difference method
a_2_bc = (dt*(Lamda_pc+beta)/2.0 + R*(1.0 - V_p*dx/D_c) + K_1 *(1+V_p*dx/D_c))
c_2_bc = K_2 + dt*(Lamda_ps+alpha)/2.0

a_2_bc_l = (dt*(Lamda_pc+beta)/2.0 + R*(1.0 + V_p*dx/D_c) + K_1 *(1-V_p*dx/D_c))

BC_PC_0 = 2.*(K_1-R)*dx*Jcin/D_c
BC_PS_0 = 2.*K_2*dx*Jsin/D_s 

# creating matrices
A1 = np.diagflat([-a_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.+a_2_bc]+[1.+a_2 for i in range(x_points-2)]+[1.+a_2_bc_l]) +\
      np.diagflat([-a_1 for i in range(x_points-1)], 1)
A2 = np.diagflat([-dt*alpha/2.0 for i in range(x_points)])
B1 = np.diagflat([a_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.-a_2_bc]+[1.-a_2 for i in range(x_points-2)]+[1.-a_2_bc_l]) +\
      np.diagflat([a_1 for i in range(x_points-1)], 1)
B2 = np.diagflat([dt*alpha/2.0 for i in range(x_points)])
C1 =  np.diagflat([-c_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.+c_2_bc]+[1.+c_2 for i in range(x_points-2)]+[1.+c_2_bc]) +\
      np.diagflat([-c_1 for i in range(x_points-1)], 1)
C2 = np.diagflat([-dt*beta/2.0 for i in range(x_points)])
D1 =  np.diagflat([c_1 for i in range(x_points-1)], -1) +\
      np.diagflat([1.-c_2_bc]+[1.-c_2 for i in range(x_points-2)]+[1.-c_2_bc]) +\
      np.diagflat([c_1 for i in range(x_points-1)], 1)
D2 = np.diagflat([dt*beta/2.0 for i in range(x_points)])

# """

"""
# BC center difference method

a_2_bc_0 = (dt*(Lamda_pc+beta)/2.0 - 2.*R* V_p*dx/D_c + 2.*K_1 *(1+V_p*dx/D_c))
a_2_bc_l = (dt*(Lamda_pc+beta)/2.0 + 2.*R*(V_p*dx/D_c) + 2*K_1 *(1-V_p*dx/D_c))
c_2_bc_0 = 2.*K_2 + dt*(Lamda_ps+alpha)/2.0
c_2_bc_l = 2.*K_2

BC_PC_0 = 4.*(K_1-R)*dx*Jcin/D_c
BC_PS_0 = 4.*K_2*dx*Jsin/D_s 

# creating matrices
A1 = np.diagflat([-a_1 for i in range(x_points-2)]+[-2.*a_1], -1) +\
      np.diagflat([1.+a_2_bc_0]+[1.+a_2 for i in range(x_points-2)]+[1.+a_2_bc_l]) +\
      np.diagflat([-2.*a_1]+[-a_1 for i in range(x_points-2)], 1)
A2 = np.diagflat([-dt*alpha/2.0 for i in range(x_points)])
B1 = np.diagflat([a_1 for i in range(x_points-2)]+[2.*a_1], -1) +\
      np.diagflat([1.-a_2_bc_0]+[1.-a_2 for i in range(x_points-2)]+[1.-a_2_bc_l]) +\
      np.diagflat([2.*a_1]+[a_1 for i in range(x_points-2)], 1)
B2 = np.diagflat([dt*alpha/2.0 for i in range(x_points)])
C1 =  np.diagflat([-c_1 for i in range(x_points-2)]+[-2.*c_1], -1) +\
      np.diagflat([1.+c_2 for i in range(x_points)]) +\
      np.diagflat([-2.*c_1]+[-c_1 for i in range(x_points-2)], 1)
C2 = np.diagflat([-dt*beta/2.0 for i in range(x_points)])
D1 =  np.diagflat([c_1 for i in range(x_points-2)] + [2.*c_1], -1) +\
      np.diagflat([1.-c_2 for i in range(x_points)]) +\
      np.diagflat([2.*c_1]+[c_1 for i in range(x_points-2)], 1)
D2 = np.diagflat([dt*beta/2.0 for i in range(x_points)])
"""



A = np.hstack((A1,A2))
# print(A)

B = np.hstack((B1,B2))

C = np.hstack((C2,C1))

D = np.hstack((D2,D1))

A_sp = sparse.csr_matrix(A)
B_sp = sparse.csr_matrix(B)
C_sp = sparse.csr_matrix(C)
D_sp = sparse.csr_matrix(D)
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
ps = np.hstack((P_c_now,P_s_now))
ps_init = ps
# breakpoint()
BC_arr= np.zeros(ps.shape)
# BC_PS_0 = np.zeros(ps.shape)
BC_arr[0] =  BC_PC_0 #influx in pc
BC_arr[len(P_s_now)] = BC_PS_0 #influx in ps
X1 = np.vstack((A,C))
Y1 = np.vstack((B,D))
X1_sp = sparse.csr_matrix(X1)
Y1_sp = sparse.csr_matrix(Y1)
breakpoint()
ps_new = ps
delta_pspine = []
def BC(ya,yb):
    return np.array([ya[1] + Jsin/D_s, yb[1],  ya[3] - V_p*ya[2]/D_c + Jcin/D_c,yb[3]- V_p*yb[2]/D_c])

def CheckBC(p_s,p_c):
    dps_a = (p_s[1]-p_s[0])/dx
    dpc_a = (p_c[1]-p_c[0])/dx
    pc_a = p_c[0]
    ps_a = p_s[0]
    ya = [ps_a,dps_a,pc_a,dpc_a]
    dps_b = (p_s[-1]-p_s[-2])/dx
    dpc_b = (p_c[-1]-p_c[-2])/dx
    pc_b = p_c[-1]
    ps_b = p_s[-1]
    yb = [ps_b,dps_b,pc_b,dpc_b]
    print("ya, yb = ",ya,yb)
    bcs = BC(ya,yb)
    return bcs
print("BC initial = ",CheckBC(P_s_init,P_c_init))

# for ti in range(1,t_points):
#     # if (ti*dt/T)*100 % t_job_percent == 0:
#     #     # breakpoint()
#     #     print(ti)
#     #     saveoutput(op_dir,date_time,P_C,P_S,P_SPINE,int((ti*dt/T)*100))
    
#     ps_new = spsolve(X1, (Y1.dot(ps) + BC_arr) )
    
#     P_c_now,P_s_now = np.split(ps_new,2)
#     delta_pspine = Spine_rk(P_s_now,P_spine_now,dt)
#     # breakpoint()
#     # print(delta_pspine[:10])
   
#     P_spine_now = P_spine_now + delta_pspine
#     P_s_now -= delta_pspine
#     ps = np.hstack((P_c_now,P_s_now))
#     residual = ps_new - ps
#     # print(residual.sum(),residual.max(),delta_pspine.sum())
#     print("Boundary conditions = ",CheckBC(P_s_now,P_c_now))
#     # res_ps = ps_new/ps_init
#     # print(res_ps[0:50],res_ps[P_s_now.shape[0]:P_s_now.shape[0]+50],delta_pspine[0:50])
#     # ps = ps_new
#     P_s_now_list = P_s_now.tolist()
#     P_c_now_list = P_c_now.tolist()
#     P_spine_now_list = P_spine_now.tolist()
#     # P_s_now = ps[200:]
#     # P_spine_now = P_spine_now
    
#     P_C.append(P_c_now_list)
#     P_S.append(P_s_now_list)
#     P_SPINE.append(P_spine_now)


# print("date and time:",date_time,"creating animations and saving")	
# saveoutput(op_dir,date_time,P_C,P_S,P_SPINE,100)
# # breakpoint()
# P_C = np.asarray(P_C)
# P_S = np.asarray(P_S)
# # P_SPINE = np.asarray(P_SPINE)

# # CallWithRunId(date_time)

