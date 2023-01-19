#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
import numpy as np
from scipy.integrate import odeint, solve_ivp
def Endocytosis(ps,alpha):
    return alpha*ps;

def Exocytosis(pc,beta):
    return pc*beta;

def Degrdation(p,lamda_p):
    return p*lamda_p

def SpineExchange(ps,pspine,eta,omega,gamma):
    return eta*ps*(omega - pspine) - gamma* pspine;


def Spine_rk(p_s, p_spine,eta,omega,gamma):
    f_vec = lambda p_s,p_spine: eta*p_s*(omega - p_spine) - gamma* p_spine
    # print(f_vec)
    k1 = f_vec(p_s, p_spine)
    k2 = f_vec(p_s + k1*dt/2.0, p_spine + k1*dt/2.0)
    k3 = f_vec(p_s + k2*dt/2.0, p_spine + k2*dt/2.0)
    k4 = f_vec(p_s + dt*k3, p_spine + k3*dt)
    
    return (1/6)* (k1 + 2.0*k2 + 2.0*k3 + k4)
def AMPATimeDynamics(t,y,Dc,Ds,Vp,lamdac,lamdas,alpha,beta,eta,omega,gamma,Jcin,Jsin,dx):
     """
    Differential equations for the 1-D coupled diffusion-advection-reaction-trapping equations.

    The ODEs are derived using the method of lines.
    """
    # The vectors pc, ps and pspineare interleaved in y.  We define
    # views of pc,ps and pspine by slicing y.
     print("time = ",t,)
     pc = y[::3]
     ps = y[1::3]
     pspine = y[2::3]
    
    # dydt is the return value of this function.
     dydt = np.empty_like(y)
    
    # Just like pc,ps and pspine are views of the interleaved vectors
    # in y, dpcdt and dpsdt,dpspinedt are views of the interleaved output
    # vectors in dydt.
     dpcdt = dydt[::3]
     dpsdt = dydt[1::3]
     dpspinedt = dydt[2::3]
     # print(y)
     # print(dpcdt.shape,dpsdt.shape,dpspinedt.shape,(Vp*(pc[1:]-pc[0:-1])/dx).shape)
     # breakpoint()
     dpcdt[0] = Endocytosis(ps[0], alpha) - Exocytosis(pc[0], beta) - Degrdation(pc[0], lamdac) \
         + 2.*Jcin*(1/dx - Vp/Dc) + Dc *(- 2.*(1+dx*Vp/Dc)*pc[0] + 2.*pc[1])/dx**2 - (Vp/dx)*((1+2*dx*Vp/Dc)*pc[0] -pc[1])
     dpcdt[1:-1] = Endocytosis(ps[1:-1], alpha) - Exocytosis(pc[1:-1], beta) - Degrdation(pc[1:-1], lamdac) \
         + Dc*np.diff(pc,2) - Vp*np.diff(pc,1)[1:]
     dpcdt[-1] = Endocytosis(ps[-1], alpha) - Exocytosis(pc[-1], beta) - Degrdation(pc[-1], lamdac) \
          + Dc*(2.*pc[-2] - 2.*(1-dx*Vp/Dc) *pc[-1])/dx**2 - (Vp/dx)*(pc[-1] - pc[-2])
          
     dpsdt[0] = Exocytosis(pc[0], beta) - Endocytosis(ps[0], alpha) - Degrdation(ps[0], lamdas)\
        + 2.*Jsin/dx + Ds*(2.*ps[1] - 2.*ps[0])/dx**2 - SpineExchange(ps[0],pspine[0],eta,omega,gamma)
     dpsdt[1:-1] = Exocytosis(pc[1:-1], beta) - Endocytosis(ps[1:-1], alpha) - Degrdation(ps[1:-1], lamdas)\
         + Ds*np.diff(ps,2) - SpineExchange(ps[1:-1],pspine[1:-1],eta,omega,gamma)
     dpsdt[-1] = Exocytosis(pc[-1], beta) - Endocytosis(ps[-1], alpha) - Degrdation(ps[-1], lamdas) \
          + Ds*(- 2.*ps[-1] + 2.*ps[-2])/dx**2 - SpineExchange(ps[-1],pspine[-1],eta,omega,gamma)
     dpspinedt[0] = SpineExchange(ps[0],pspine[0],eta,omega,gamma)
     dpspinedt[1:-1] = SpineExchange(ps[1:-1],pspine[1:-1],eta,omega,gamma) #Spine_rk(ps,pspine,eta,omega,gamma)
     dpspinedt[-1] = SpineExchange(ps[-1],pspine[-1],eta,omega,gamma)
     # dpsdt -= dpspinedt
     # print(dpcdt,dpsdt,dpspinedt,"*"*20,"\n")
     return dydt
     


from AMPA_model import RunSim5, RunModelWithFile
from datetime import datetime
import matplotlib.pyplot as plt
P_s_init,P_c_init,P_spine_init,SP_model1 = RunModelWithFile("./ModelParams.json")
y0 = np.vstack((np.vstack((P_c_init,P_s_init)),P_spine_init))
L = 500.0
dx = SP_model1.dx
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
y_init_orig = y0.T.flatten()#np.zeros(y0.T.flatten().shape)
dt = 0.002
t_max = 2
ps= y_init_orig[::3]
pc= y_init_orig[1::3]
pspine = y_init_orig[2::3]
x_grid = np.arange(0,L,dx)
y_init = y_init_orig
def BC(ya,yb):
    return np.array([ya[1] + Jsin/D_s, yb[1],  ya[3] - V_p*ya[2]/D_c + Jcin/D_c,yb[3]- V_p*yb[2]/D_c])

def CheckBC(p_s,p_c):
    dps_a = np.diff(p_s,1)[0]
    dpc_a = np.diff(p_c,1)[0]
    pc_a = p_c[0]
    ps_a = p_s[0]
    ya = [ps_a,dps_a,pc_a,dpc_a]
    dps_b = np.diff(p_s,1)[-1]
    dpc_b = np.diff(p_c,1)[-1]
    pc_b = p_c[-1]
    ps_b = p_s[-1]
    yb = [ps_b,dps_b,pc_b,dpc_b]
    print("ya, yb = ",ya,yb)
    bcs = BC(ya,yb)
    return bcs
print("BC initial = ",CheckBC(P_s_init,P_c_init))
print("Integal condition = ",SP_model1.IntegralBC(P_s_init, P_c_init))
# breakpoint()
t_steps = int(t_max/dt)
# for i in range(0,t_steps):
#     t_step = [i*dt,(i+1)*dt]
#     # breakpoint()
#     try:
#         soln = solve_ivp(AMPATimeDynamics,t_step,y_init,args=(D_c,D_s,V_p,Lamda_pc,Lamda_ps,alpha,beta,eta,omega,gamma,Jcin,Jsin,dx), dense_output=True,vectorize=True,lband=2,uband=2,rtol=1,atol=1,hmax=dt)
#     except:
#         breakpoint()
#         print("error occured")
#     y_init = soln.y[:-1]
#     if np.isnan(y_init).any():
#         print("first nana encountered")
#         breakpoint()
# soln = y_init.y
t = np.arange(0,t_max,dt)
soln = solve_ivp(AMPATimeDynamics,[0,t_max],y_init,args=(D_c,D_s,V_p,Lamda_pc,Lamda_ps,alpha,beta,eta,omega,gamma,Jcin,Jsin,dx), dense_output=True,vectorize=True,lband=2,uband=2,rtol=1e-5,atol=1e-11,max_step=dt)
breakpoint()
pc_final = soln.y[::3,-1]
ps_final = soln.y[1::3,-1]
pspine_final = soln.y[2::3,-1]
print("BC initial = ",CheckBC(ps_final,pc_final))
print("Integal condition = ",SP_model1.IntegralBC(ps_final, pc_final))
plt.plot(x_grid,pc_final,'r:',label = r'$p_{c\infty+100}$')
plt.plot(x_grid,ps_final,'g:',label = r'$p_{s\infty+100}$')
plt.plot(x_grid,pspine_final,'b:',label = r'$p_{spine\infty+100}$')
plt.plot(x_grid,P_c_init,'r--',label = r'$p_{c\infty}$')
plt.plot(x_grid,P_s_init,'g--',label = r'$p_{s\infty}$')
plt.plot(x_grid,P_spine_init,'b--',label = r'$p_{spine\infty}$')
plt.legend()
plt.show()


breakpoint()