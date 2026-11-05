#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
import numpy as np
from TemporalIntegration import *
from datetime import datetime

pc = y_init_orig[::3]
ps = y_init_orig[1::3]
pspine = y_init_orig[2::3]
x_grid = np.arange(0, L, dx)
y_init = y_init_orig

def _1gaussian(x,cen1,sigma1):

    return (np.exp((-((x-cen1)/sigma1)**2)))
# def getNet_interval(locations,span):
def Getchange_profile(x_grid,locns,sigma,factor):
    change_profile = np.ones(x_grid.shape)
    for l1 in locns:
        print(l1,sigma,factor)
        change_profile += (factor*_1gaussian(x_grid,l1,sigma))
    # plt.plot(x_grid,change_profile)
    # plt.show()
        # breakpoint()
    return change_profile
def ParamChangeGauss(x_grid,param, locns, sigma, factor,uod):
    bf = factor
    b_sig= sigma
    beta= param
    beta_final = beta
#     for each parameter, if the change factor is not = 1, we put a gaussian change at
#     each location in locns array with mu = location, sigma = sigma and amp = factor
    if not bf ==1:
        # beta_final = beta
        b_change_profile = Getchange_profile(x_grid,locns,b_sig,bf)
        if uod == 1:
            beta_final  = beta*b_change_profile
        elif uod == -1:
            beta_final = beta / b_change_profile
        else:
            return beta

    return beta_final

def Singledxchange(pa,locns,factor,uod,dx):
    p = pa.copy()
    for l1 in locns:
        p[int(l1/dx)]  *= ((factor)**uod)
    return p
def PlasticityExperimentGauss(x_grid,params,param_names, locns, sigmas, factors,up_or_down,step_num,op_dir,dat_time,dx):
    new_pa = params.copy()
    fig,ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    plt.yscale("log")
    for pdx,p in enumerate(params):
        if sigmas[pdx] == dx:
            # print("modifying ",param_names[pdx],((factors[pdx])**up_or_down[pdx]))
            # for l1 in locns:
            new_pa[pdx] = Singledxchange(params[pdx],locns,factors[pdx],up_or_down[pdx],dx)
            # for l1 in locns:
            #     print( new_pa[pdx][int(l1/dx)], new_pa[pdx][int(l1/dx)-1],params[pdx][int(l1/dx)], params[pdx][int(l1/dx)-1])
        else:
            new_pa[pdx] = ParamChangeGauss(x_grid,params[pdx],locns,sigmas[pdx],factors[pdx],up_or_down[pdx])
        if not factors[pdx] == 1:
            ax.plot(x_grid, new_pa[pdx] / params[pdx], label=param_names[pdx])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(r"Distance from soma ($\mu$m)")
    ax.set_ylabel(r"$Log_{10}$[Fold change]")
    ax.tick_params(axis='both', which='major', labelsize=18)
    # plt.xlabel("Simulation time in mins")
    # plt.ylabel("Fold change")
    plt.legend(frameon=False, fontsize=18)

    # SaveFigures("{0}/fig_protocol_{1}_step_{2}".format(op_dir, dat_time, step_num), dpi=300)
    plt.show()
    return new_pa


def CNIH2dependentGluA2Plasticity(ds,dc,vp,alpha,beta_0,beta,eta,gamma,lo,c2_fname):
    """
    Model of CNIH2 local translation dependent plasticity
    CNIH2 local translation drives the change in exocytosis
    Change is spine size lead to longer dwell time in synapses hence slower unbinding rate for P_spine
    """
    sim_type = "G2_CNIH2_dep_cLTP"
    sim_time = 0
    time_steps = []

    b_factors = []
    a_factors = []
    g_factors = []
    dc_factors = []
    vp_factors = []
    ds_factors = []
    eta_factors = []

    l_scales0 = [dx,dx,dx,dx,dx,dx]
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/Time-dependent/GluA2/CNIH2-dep/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    # breakpoint()
    """
    Step 0
    """
    # f0 = 1  # increase by a factor of 1
    # a0 = 1
    g0 = 1
    # dc0 = 1
    # vp0 = 1
    # ds0 = 1
    eta0 = 1.5 

    t_step0 = 5*60*60  # t secs
    # b_factors.append(f0)
    # a_factors.append(a0)
    g_factors.append(g0)
    eta_factors.append(eta0)
    # dc_factors.append(dc0)
    # vp_factors.append(vp0)
    # ds_factors.append(ds0)
    time_steps.append(t_step0)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)
    eta1 = PlasticityExperimentGauss(x_grid,[eta],["eta"], lo, [dx], [eta0],[1], 1,"",date_time,dx)[0]
    # breakpoint()
    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha , beta_0,beta, eta1, omega, gamma*g0,
                    Jcin, Jsin, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = AdaptiveDynamicsSimRun(model_params, t_range, t_eval, y_init, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)
#    breakpoint()
    
    meta_dir = os.getcwd() + "/Time-dependent/GluA2/CNIH2-dep/meta_" + date_time
    os.makedirs(meta_dir,exist_ok=True)
    print("date and time:", date_time)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    savesimsettings(1,sim_type, time_steps,"{0}/protocol_{1}.json".format(meta_dir, date_time),beta_factors=b_factors, eta_factors = eta_factors,locations = lo.tolist(),
                    gamma_factors=g_factors, dc_factors=dc_factors,vp_factors=vp_factors,ds_factors=ds_factors,alpha_factors=a_factors,cnih2_file = c2_fname)

def CNIH2IndependentGluA2Plasticity(ds,dc,vp,alpha,beta_0,eta,gamma,lo,global_ltp = True):
    """
    Model of GluA2 plasticity response without inducing CNIH2 local translation
    step 0 is as always, running the simulation without any parameter changes for 30 secs
    step1 is : parameters are changed and the simulation is run for time t_step1
    step 2 is same as step 1 with different changes 
    finally the results are saved
    """
    sim_type = "G2_CNIH2_Indep_cLTP"
    sim_time = 0
    time_steps = [0]
    e_factor = []
    b_factor = []
    ds_factor = []
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    """
    Step 0
    """
    e0 = 1
    b0 = 1
    ds0 = 1
    t_step0 = 30  # t secs
    e_factor.append(e0)
    b_factor.append(b0)
    ds_factor.append(ds0)
    time_steps.append(t_step0)
    model_params = [dc, ds*ds0, vp, Lamda_pc, Lamda_ps, alpha , beta_0*b0, eta*e0, omega, gamma,
                    Jcin, Jsin, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, y_init, max_step= dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)

    """
    step1
    """
    ds1 = 1
    e1_f = 2 
    b1_f = 10
    if global_ltp:
        
        b1 = 10 
        model_params = [dc, ds*ds1, vp, Lamda_pc, Lamda_ps, alpha, beta_0*b1, eta * e1, omega, gamma,
                        Jcin, Jsin, dx]
    else:
        b1,e1 = PlasticityExperimentGauss(x_grid,[beta_0,eta],["beta","eta"], lo, [dx,dx], [b1_f,e1_f],[1,1], 1,"",date_time,dx)
        model_params = [dc, ds*ds1, vp, Lamda_pc, Lamda_ps, alpha, b1,  e1, omega, gamma,
                        Jcin, Jsin, dx]
    # breakpoint()
  
    t_step1 = 2*60*60  # t secs
    e_factor.append(e1_f)
    b_factor.append(b1_f)
    ds_factor.append(ds1)
    time_steps.append(t_step1)
    t_range = [0, t_step1]
    t_eval = np.arange(0, t_step1, dt)
    y_step1 = data_mat[:,-1]
    soln1 = DynamicSimRun(model_params, t_range, t_eval, y_step1, max_step=100 * dt, method='RK45')
    data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    total_tps = np.concatenate((total_tps,(soln1.t + sim_time)))
    sim_time += t_step1
    print("Step 1 finished at simulation time  = ", sim_time)
    """
    step2
    """
    e2 = e1_f
    b2 =  1
    ds2 = 1 
    t_step2 = 5*60*60 - t_step1  # t secs
    e_factor.append(e2)
    b_factor.append(b2)
    ds_factor.append(ds2)
    time_steps.append(t_step2)
    model_params = [dc, ds*ds2, vp, Lamda_pc, Lamda_ps, alpha, beta_0*b2, eta * e2, omega, gamma,
                    Jcin, Jsin, dx]
    t_range = [0, t_step2]
    t_eval = np.arange(0, t_step2, dt)
    y_step2 = data_mat[:,-1]
    soln2 = DynamicSimRun(model_params, t_range, t_eval, y_step2, max_step=100 * dt, method='RK45')
    data_mat = np.concatenate((data_mat, soln2.y), axis=1)
    total_tps = np.concatenate((total_tps,(soln2.t + sim_time)))
    sim_time += t_step2
    print("Step 2 finished at simulation time  = ", sim_time)


    """
    saving everything
    """
    
    op_dir = os.getcwd() + "/Time-dependent/GluA2/CNIH2-Indep/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    savesimsettings(3, sim_type, time_steps,"{0}/protocol_{1}.json".format(op_dir, date_time),ds_factors = ds_factor,beta_factors=b_factor,eta_factors=e_factor)
    print("Simulation finished without any error")

CNIH2_dep = True 

beta_0 = beta_orig
gamma_arr = gamma_orig*np.ones(P_c_init.shape)
dc_arr = D_c_orig*np.ones(P_c_init.shape)
vp_arr = V_p_orig*np.ones(P_c_init.shape)
ds_arr = D_s_orig*np.ones(P_c_init.shape)
alpha_arr = alpha_orig*np.ones(P_c_init.shape)
eta_arr = eta_orig*np.ones(P_c_init.shape)
n = 10
np.random.seed(2026)
location = np.random.randint(10,501,size =n)
# x_span_dx = int(x_span )
# prolonged local: 11_08_2024_16_13_52,04_14_2026_20_45_18
# brief local: 11_11_2024_10_24_43 
# somatic : 11_11_2024_11_56_23
if CNIH2_dep:
    cwd = os.getcwd()
    CNIH2_rf = os.path.join(cwd,"CNIH2-translation")
    print(CNIH2_rf)
    date_n_time = "11_08_2024_16_13_52"
    data_file = "CNIH2_{}_100_percent.npy".format(date_n_time)
    folder_to_Read = os.path.join(CNIH2_rf, date_n_time)
    beta_mat = np.load(os.path.join(folder_to_Read, data_file))
    CNIH2dependentGluA2Plasticity(ds_arr,dc_arr,vp_arr,alpha_arr,beta_0,beta_mat,eta_arr,gamma_arr, location,date_n_time)
else:
    beta_arr = beta_orig*np.ones(P_c_init.shape)
    CNIH2IndependentGluA2Plasticity(ds_arr, dc_arr, vp_arr, alpha_arr, beta_arr, eta_arr, gamma_arr, location,global_ltp = False)

