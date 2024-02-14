#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
from TemporalIntegration import *
from datetime import datetime

pc = y_init_orig[::3]
ps = y_init_orig[1::3]
pspine = y_init_orig[2::3]
x_grid = np.arange(0, L, dx)
y_init = y_init_orig


def CNIH2dependentGluA2Plasticity(ds,dc,vp,alpha,beta_0,beta,eta,gamma,lo):
    """
    Model of CNIH2 local translation dependent plasticity
    CNIH2 local translation drives the change in exocytosis
    Change is spine size lead to longer dwell time in synapses hence slower unbinding rate for P_spine
    """

    sim_time = 0
    time_steps = []

    b_factors = []
    a_factors = []
    g_factors = []
    dc_factors = []
    vp_factors = []
    ds_factors = []

    l_scales0 = [dx,dx,dx,dx,dx,dx]
    # breakpoint()
    """
    Step 0
    """
    # f0 = 1  # increase by a factor of 1
    # a0 = 1
    g0 = 1/1.5
    # dc0 = 1
    # vp0 = 1
    # ds0 = 1

    t_step0 = 5*60*60  # t secs
    # b_factors.append(f0)
    # a_factors.append(a0)
    g_factors.append(g0)
    # dc_factors.append(dc0)
    # vp_factors.append(vp0)
    # ds_factors.append(ds0)
    time_steps.append(t_step0)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)
    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha , beta_0,beta, eta, omega, gamma*g0,
                    Jcin, Jsin, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = AdaptiveDynamicsSimRun(model_params, t_range, t_eval, y_init, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)
    breakpoint()
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/Time-dependent/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    savesimsettings(1, time_steps, lo,"{0}/protocol_{1}.json".format(op_dir, date_time),beta_factors=b_factors,
                    gamma_factors=g_factors, dc_factors=dc_factors,vp_factors=vp_factors,ds_factors=ds_factors,alpha_factors=a_factors)

def CNIH2IndependentGluA2Plasticity(ds,dc,vp,alpha,beta_0,eta,gamma,lo):
    """
    Model of CNIH2 local translation dependent plasticity
    CNIH2 local translation drives the change in exocytosis
    Change is spine size lead to longer dwell time in synapses hence slower unbinding rate for P_spine
    """

    sim_time = 0
    time_steps = []

    b_factors = []
    a_factors = []
    g_factors = []
    dc_factors = []
    vp_factors = []
    ds_factors = []

    l_scales0 = [dx,dx,dx,dx,dx,dx]
    # breakpoint()
    """
    Step 0
    """
    # f0 = 1  # increase by a factor of 1
    # a0 = 1
    g0 = 1/1.5
    # dc0 = 1
    # vp0 = 1
    # ds0 = 1

    t_step0 = 1*60*60  # t secs
    # b_factors.append(f0)
    # a_factors.append(a0)
    g_factors.append(g0)
    # dc_factors.append(dc0)
    # vp_factors.append(vp0)
    # ds_factors.append(ds0)
    time_steps.append(t_step0)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)
    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha , beta_0, eta, omega, gamma*g0,
                    Jcin, Jsin, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, y_init, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)
    breakpoint()
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/Time-dependent/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    savesimsettings(1, time_steps, lo,"{0}/protocol_{1}.json".format(op_dir, date_time),beta_factors=b_factors,
                    gamma_factors=g_factors, dc_factors=dc_factors,vp_factors=vp_factors,ds_factors=ds_factors,alpha_factors=a_factors)
    print("Simulation finished without any error")

CNIH2_dep = False

beta_0 = beta_orig
gamma_arr = gamma_orig*np.ones(P_c_init.shape)
dc_arr = D_c_orig*np.ones(P_c_init.shape)
vp_arr = V_p_orig*np.ones(P_c_init.shape)
ds_arr = D_s_orig*np.ones(P_c_init.shape)
alpha_arr = alpha_orig*np.ones(P_c_init.shape)
eta_arr = eta_orig*np.ones(P_c_init.shape)
loc = [250]
location = [int(l) for l in loc]
# x_span = 3
# x_span_dx = int(x_span )
if CNIH2_dep:
    CNIH2_rf = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Anne-Sophie/CNIH2-translation/"
    date_n_time = "10_17_2023_16_24_23"
    data_file = "CNIH2_{}_100_percent.npy".format(date_n_time)
    folder_to_Read = os.path.join(CNIH2_rf, date_n_time)
    beta_mat = np.load(os.path.join(folder_to_Read, data_file))
    CNIH2dependentGluA2Plasticity(ds_arr,dc_arr,vp_arr,alpha_arr,beta_0,beta_mat,eta_arr,gamma_arr, location)
else:
    beta_arr = beta_orig*np.ones(P_c_init.shape)
    CNIH2IndependentGluA2Plasticity(ds_arr, dc_arr, vp_arr, alpha_arr, beta_arr, eta_arr, gamma_arr, location)

