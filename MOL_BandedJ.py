#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
from AMPA_model import RunSim5, RunModelWithFile
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.integrate import odeint, solve_ivp
import shutil
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
     if t%2 == 0:
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
     # print(dpcdt.shape,dpsdt.shape,dpspinedt.shape,pc.shape,ps.shape,pspine.shape)
     # breakpoint()
     dpcdt[0] = Endocytosis(ps[0], alpha[0]) - Exocytosis(pc[0], beta[0]) - Degrdation(pc[0], lamdac) \
         + 2.*Jcin*(1/dx - Vp/Dc) + Dc *(- 2.*(1+dx*Vp/Dc)*pc[0] + 2.*pc[1])/dx**2 - (Vp/dx)*((1+2*dx*Vp/Dc)*pc[0] -pc[1])
     dpcdt[1:-1] = Endocytosis(ps[1:-1], alpha[1:-1]) - Exocytosis(pc[1:-1], beta[1:-1]) - Degrdation(pc[1:-1], lamdac) \
         + Dc*np.diff(pc,2) - Vp*np.diff(pc,1)[1:]
     dpcdt[-1] = Endocytosis(ps[-1], alpha[-1]) - Exocytosis(pc[-1], beta[-1]) - Degrdation(pc[-1], lamdac) \
          + Dc*(2.*pc[-2] - 2.*(1-dx*Vp/Dc) *pc[-1])/dx**2 - (Vp/dx)*(pc[-1] - pc[-2])
          
     dpsdt[0] = Exocytosis(pc[0], beta[0]) - Endocytosis(ps[0], alpha[0]) - Degrdation(ps[0], lamdas)\
        + 2.*Jsin/dx + Ds*(2.*ps[1] - 2.*ps[0])/dx**2 - SpineExchange(ps[0],pspine[0],eta,omega,gamma)
     dpsdt[1:-1] = Exocytosis(pc[1:-1], beta[1:-1]) - Endocytosis(ps[1:-1], alpha[1:-1]) - Degrdation(ps[1:-1], lamdas)\
         + Ds*np.diff(ps,2) - SpineExchange(ps[1:-1],pspine[1:-1],eta,omega,gamma)
     dpsdt[-1] = Exocytosis(pc[-1], beta[-1]) - Endocytosis(ps[-1], alpha[-1]) - Degrdation(ps[-1], lamdas) \
          + Ds*(- 2.*ps[-1] + 2.*ps[-2])/dx**2 - SpineExchange(ps[-1],pspine[-1],eta,omega,gamma)
     dpspinedt[0] = SpineExchange(ps[0],pspine[0],eta,omega,gamma)
     dpspinedt[1:-1] = SpineExchange(ps[1:-1],pspine[1:-1],eta,omega,gamma) #Spine_rk(ps,pspine,eta,omega,gamma)
     dpspinedt[-1] = SpineExchange(ps[-1],pspine[-1],eta,omega,gamma)
     dydt
     # dpsdt -= dpspinedt
     # print(dpcdt,dpsdt,dpspinedt,"*"*20,"\n")
     return dydt
     

def saveoutput(op_dir,date_time,soln,tps,percent,orig_params_file):
    p_c,p_s,p_sp = GetProteindata(soln)
    with open('{0}/PC_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_c)
    f.close()
    with open('{0}/PS_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_s)
    f.close()
    with open('{0}/PSPINE_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,p_sp)
    f.close()
    t_points = tps
    with open('{0}/timepoints_{1}_{2}_percent.npy'.format(op_dir,date_time,percent), 'wb') as f:
        np.save(f,t_points)
    f.close()
    shutil.copyfile(orig_params_file, '{0}/baseline_parameters_{1}.json'.format(op_dir,date_time))

def GetProteindata(soln):
    pc_t = soln[::3]
    ps_t = soln[1::3]
    p_sp_t = soln[2::3]
    return pc_t,ps_t,p_sp_t
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
y_init_orig = y0.T.flatten()#np.zeros(y0.T.flatten().shape) #np.zeros(y0.T.flatten().shape) #
dt = 0.02
t_max = 2000
times =  100
t = np.arange(0,t_max,times*dt)  #times for which to save the y values
pc= y_init_orig[::3]
ps= y_init_orig[1::3]
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
def DynamicSimRun(model_params,t_range,t_eval,y_init,method="LSODA",dense_op=True,vectorize=True,lband=2,uband=2,rtol=1e-5,atol=1e-11,max_step=100):
    soln = solve_ivp(AMPATimeDynamics,t_range,y_init,args=(model_params), \
                     method=method,dense_output=dense_op,vectorize=vectorize,lband=lband,uband=uband,rtol=rtol,atol=atol,max_step=max_step,t_eval = t_eval)
    return soln
def PlasticityExperiment(beta,lo,x_sdx,up_or_down_factor):
    beta_updated = beta.copy()
    print(up_or_down_factor)
    beta_updated[lo-x_sdx:lo+x_sdx] *= up_or_down_factor
    return beta_updated


def PattersonGluA1Plasticity(beta,lo,x_sdx):
    """
    step zero, we simulate the model in steady state for 10 secs
    step First, we step-increase the exocytosis rate (beta) to f1 folds in 50 seconds
    step Second, we step-increase the exocytosis rate (beta) to f2 folds for 10 seconds
    step last, we change back the parameters to basal level and run simulation for 30 mins
    we integrate for a total of 30 mins to see the GluA1 dynamics
    """
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/Time-dependent/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    sim_time = 0
    time_steps = [0]
    factors = [1]
    """
    Step 0
    """
    f0 = 1 #increase by a factor of 1
    t_step0 = 30  # running for 50 secs
    factors.append(f0)
    time_steps.append(t_step0)
    beta_step0 = PlasticityExperiment(beta, lo,x_sdx, f0)

    model_params = [D_c, D_s, V_p, Lamda_pc, Lamda_ps, alpha * np.ones(P_c_init.shape), beta_step0, eta, omega, gamma,
                    Jcin, Jsin, dx]
    plt.plot(x_grid, beta_step0 / beta_step0[0])
    plt.show()
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, y_init, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)
    beta_profile = np.ones(soln0.t.shape)*beta_step0[location]
    """
    Step 1
    """
    f1 = 26
    beta_step1  = PlasticityExperiment(beta, lo,x_sdx, f1)
    t_step1 = 60       #running for 50 secs
    factors.append(f1)
    time_steps.append(t_step1)
    model_params = [D_c,D_s,V_p,Lamda_pc,Lamda_ps,alpha*np.ones(P_c_init.shape),beta_step1,eta,omega,gamma,Jcin,Jsin,dx]
    plt.plot(x_grid,beta_step1/beta_step1[0])
    plt.show()
    t_range = [0,t_step1]
    t_eval = np.arange(0, t_step1, dt)
    soln1=DynamicSimRun(model_params,t_range,t_eval,y_init,max_step=100*dt,method='RK45')
    total_tps = np.concatenate((total_tps,(soln1.t+sim_time)))
    sim_time += t_step1
    data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    print("Step 1 finished at simulation time  = ", sim_time)
    saveoutput(op_dir, date_time, data_mat, total_tps, 1, "./ModelParams.json")
    beta_profile = np.concatenate((beta_profile ,np.ones(soln1.t.shape)*beta_step1[location]))
    # """
    # Step 2
    # """
    # f2 = 6
    # beta_step2 = lasticityExperiment(beta, lo,x_sdx, f2)
    # t_step2 = 10        #running for 10 secs
    # factors.append(f2)
    # time_steps.append(t_step2)
    # # breakpoint()
    # """
    # chaning y_init to the last time_step value of step 1
    # """
    # new_y_init = soln1.y[:,-1]
    # model_params2 = [D_c, D_s, V_p, Lamda_pc, Lamda_ps, alpha * np.ones(P_c_init.shape), beta_step2, eta, omega, gamma,
    #                 Jcin, Jsin, dx]
    # plt.plot(x_grid, beta_step2/beta_step2[0])
    # plt.show()
    # t_eval2 = np.arange(0, t_step2, dt)
    # t_range2 = [0, t_step2]
    # soln2 = DynamicSimRun(model_params2,t_range2,t_eval2,new_y_init,max_step=100*dt,method='RK45')
    #
    # data_mat = np.concatenate((data_mat, soln2.y), axis=1)
    # total_tps = np.concatenate((total_tps, (soln2.t + sim_time)))
    # saveoutput(op_dir, date_time, data_mat, total_tps, 2, "./ModelParams.json")
    # sim_time += t_step2
    # print("Step 2 finished at simulation time  = ", sim_time)
    # beta_profile = np.concatenate((beta_profile, np.ones(soln2.t.shape) * beta_step2[location]))
    """
    Step last
    """
    finf = 1
    beta_steplast  = PlasticityExperiment(beta, lo,x_sdx, finf)
    t_steplast = 10  # running for 10 secs
    factors.append(finf)
    time_steps.append(t_steplast)
    # breakpoint()
    """
    chaning y_init to the last time_step value of step 2
    """
    new_y_init = data_mat[:, -1]
    model_params2 = [D_c, D_s, V_p, Lamda_pc, Lamda_ps, alpha * np.ones(P_c_init.shape), beta_steplast, eta, omega, gamma,
                     Jcin, Jsin, dx]
    plt.plot(x_grid, beta_steplast / beta_steplast[0])
    plt.show()
    t_evallast= np.arange(0, t_steplast, dt)
    t_rangelast = [0, t_steplast]
    soln_last = DynamicSimRun(model_params2, t_rangelast, t_evallast, new_y_init, max_step=100 * dt, method='RK45')
    print("Last finished at simulation time  = ", sim_time)
    # saveoutput(op_dir, date_time, data_mat, total_tps, n, "./ModelParams.json")
    beta_profile = np.concatenate((beta_profile, np.ones(soln_last.t.shape) * beta_steplast[location]))
    """
    Concatinating the results of all steps
    """
    total_tps = np.concatenate((total_tps, (soln_last.t + sim_time)))
    sim_time += t_steplast
    print("total simulation time = ",sim_time)
    data_mat = np.concatenate((data_mat, soln_last.y),axis=1)
    saveoutput(op_dir, date_time, data_mat,total_tps, 100,"./ModelParams.json")
    plt.plot(total_tps,beta_profile)
    breakpoint()
    savesimsettings(3,time_steps,factors,[loc],x_span,"{0}/protocol_{1}".format(op_dir,date_time))
    plt.show()

def savesimsettings(num_steps,time_steps,factors,locations,x_span,protocol_file):
    assert num_steps == len(time_steps) - 1
    assert num_steps == len(factors) - 1
    with open(protocol_file, 'a') as the_file:
        the_file.write('#Step = {}\n'.format(num_steps))
        the_file.write("time points = ")
        for t in time_steps:
            the_file.write(str(t))
            the_file.write("\t")
        the_file.write("\n")
        the_file.write("factors = ")
        for f in factors:
            the_file.write(str(f))
            the_file.write("\t")
        the_file.write("\n")
        the_file.write("location = ")
        for l in locations:
            the_file.write(str(l))
            the_file.write("\t")
        the_file.write("\n")
        the_file.write("x span = {}\n".format(x_span))

        the_file.write("parameter chaned = beta")
        the_file.close()
beta_array = beta*np.ones(P_c_init.shape)
loc = 100
location = int(loc/dx)
x_span = 3
x_span_dx = int(x_span/dx)
PattersonGluA1Plasticity(beta_array,location,x_span_dx)
