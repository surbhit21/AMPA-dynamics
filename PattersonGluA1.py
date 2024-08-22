#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
import matplotlib.pyplot as plt
import numpy as np
from AMPA_model import RunSimGluA1
# RunSimGluA1(0.24,1e-5,0.1,0.1)
from TemporalIntegration import *
from datetime import datetime


# t_max = 2000
# times = 100
# t = np.arange(0, t_max, times * dt)  # times for which to save the y values
pc = y_init_orig[::3]
ps = y_init_orig[1::3]
pspine = y_init_orig[2::3]
x_grid = np.arange(0, L, dx)
y_init = y_init_orig

breakpoint()




# print("BC initial = ", CheckBC(P_s_init, P_c_init))
# print("Integal condition = ", SP_model1.IntegralBC(P_s_init, P_c_init))


# breakpoint()

def PlasticityExperiment(beta,alpha,gamma,dc,vp, ds, locns, x_sdx, up_or_down_factor,a_factor,g_factor,dc_factor,vp_factor,ds_factor,step_num,op_dir,dt):

    # change_profile = np.arange(lo-x_sdx,lo+x_sdx,1)
    print(up_or_down_factor,g_factor,dc_factor,vp_factor)
    beta_final = beta.copy()
    alpha_final = alpha.copy()
    gamma_final = gamma.copy()
    dc_final = dc.copy()
    vp_final = vp.copy()
    ds_final = ds.copy()
    for lo in locns:
        beta_updated = beta.copy()
        alpha_updated = alpha.copy()
        gamma_updated = gamma.copy()
        dc_updated = dc.copy()
        vp_updated = vp.copy()
        ds_updated = ds.copy()
        if not up_or_down_factor == 1:
            beta_updated[lo-x_sdx:lo+x_sdx] *= (up_or_down_factor)
            beta_final[lo-x_sdx:lo+x_sdx] += beta_updated[lo-x_sdx:lo+x_sdx]
        if not a_factor == 1:
            alpha_updated[lo - x_sdx:lo + x_sdx] *= (a_factor)
            alpha_final[lo - x_sdx:lo + x_sdx] += alpha_updated[lo - x_sdx:lo + x_sdx]
        if not g_factor == 1:
            gamma_updated[lo] *= g_factor
            gamma_final[lo] = gamma_updated[lo]
        # breakpoint()
        if not dc_factor == 1:
            dc_updated[lo] *= (dc_factor)
            dc_final[lo] = dc_updated[lo]
        if not vp_factor == 1:
            vp_updated[lo] *= (vp_factor)
            vp_final[lo] = vp_updated[lo]
        if not ds_factor == 1:
            ds_updated[lo] *= (ds_factor)
            ds_updated[lo] = ds_updated[lo]

    # vp_updated[lo + x_sdx] *= vp_factor
    plt.plot(x_grid, beta_final / beta, label=r"$\beta$")
    plt.plot(x_grid, alpha_final / alpha, label=r"$\alpha$")
    plt.plot(x_grid, gamma_final / gamma, label=r"$\gamma$")
    # plt.plot(x_grid, dc_final / dc, label=r"$D_c$")
    # plt.plot(x_grid, vp_final / vp, label=r"$V_p$")
    # plt.plot(x_grid, ds_final / ds, label=r"$D_s$")
    plt.title("at time step {}".format(step_num))
    plt.legend()
    SaveFigures("{0}/fig_protocol_{1}_step_{2}".format(op_dir, dt,step_num),dpi=300)
    plt.show()
    return beta_final,alpha_final, gamma_final, dc_final, vp_final, ds_final



def patterson2010GluA1Plasticity(ds,dc,vp,alpha,beta,eta,gamma, lo):
    """
    step zero, we simulate the model in steady state for 10 secs
    step First, we step-increase the exocytosis rate (beta) to f1 folds in 50 seconds
    step Second, we step-increase the exocytosis rate (beta) to f2 folds for 10 seconds
    step last, we change back the parameters to basal level and run simulation for 30 mins
    we integrate for a total of 30 mins to see the GluA1 dynamics
    """
    type = "G1_cLTP"
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
    f0 = 1  # increase by a factor of 1
    a0 = 1
    g0 = 1
    dc0 = 1
    vp0 = 1
    ds0 = 1

    t_step0 = 30  # running for t_step0 secs
    b_factors.append(f0)
    a_factors.append(a0)
    g_factors.append(g0)
    dc_factors.append(dc0)
    vp_factors.append(vp0)
    ds_factors.append(ds0)
    time_steps.append(t_step0)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/Time-dependent/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    beta_step0, alpha_Step0, gamma_step0, dc_step0, vp_step0, ds_step0 = PlasticityExperimentGauss(x_grid, [beta_array, alpha_arr, gamma_arr, dc_arr, vp_arr, ds_arr],
                              [r"$\beta$", r"$\alpha$", r"$\gamma$", r"$D_c$", r"$V_p$", r"$D_s$"], location,
                              l_scales0, [f0, a0, g0, dc0, vp0, ds0], [1, 1, 1, 1, 1, 1], 0, op_dir,date_time)
    model_params = [dc_step0 , ds_step0, vp_step0, Lamda_pc, Lamda_ps, alpha_Step0 , beta_step0, eta, omega, gamma_step0,
                    Jcin, Jsin, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, y_init, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)
    """
    Step 1
    """
    f1 = 2600
    a1 = 1
    g1 = 3
    dc1 = 1
    vp1  = 1
    ds1 = 1
    new_y_init = data_mat[:, -1]  # Bleach(data_mat[:, -1],lo,x_sdx)
    l_scales1 = [6,dx,dx,dx,dx,dx]
    # l_scales =
    # beta_step1,alpha_step1,gamma_step1,dc_step1, vp_step1, ds_step1 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f1,a1,g1,dc1,vp1,ds1,1,op_dir,date_time)
    beta_step1,alpha_step1,gamma_step1,dc_step1, vp_step1, ds_step1 = PlasticityExperimentGauss(x_grid, [beta_array, alpha_arr, gamma_arr, dc_arr, vp_arr, ds_arr],
                              [r"$\beta$", r"$\alpha$", r"$\gamma$", r"$D_c$", r"$V_p$", r"$D_s$"], location,
                              l_scales1, [f1, a1, g1, dc1, vp1, ds1], [1, 1, -1, 1, 1, 1], 1, op_dir, date_time)
    plt.plot(x_grid,beta_step1/beta_step0)
    plt.show()
    t_step1 = 1*60  # running for 60 secs
    b_factors.append(f1)
    a_factors.append(a1)
    g_factors.append(g1)
    dc_factors.append(dc1)
    vp_factors.append(vp1)
    ds_factors.append(ds1)
    time_steps.append(t_step1)
    model_params = [dc_step1 , ds_step1, vp_step1, Lamda_pc, Lamda_ps, alpha_step1 , beta_step1, eta, omega, gamma_step1,
                    Jcin, Jsin, dx]

    t_range = [0, t_step1]
    t_eval = np.arange(0, t_step1, dt)
    soln1 = DynamicSimRun(model_params, t_range, t_eval, new_y_init, max_step=100 * dt, method='RK45')
    total_tps = np.concatenate((total_tps, (soln1.t + sim_time)))
    sim_time += t_step1
    data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    print("Step 1 finished at simulation time  = ", sim_time)
    saveoutput(op_dir, date_time, data_mat, total_tps, 1, baseline_param_file)
    """
    Step last
    """
    finf = 1
    ainf = 1
    ginf = 2
    dcinf = 1
    vpinf = 1
    dsinf = 1
    l_scalesinf = l_scales0
    # beta_steplast,alpha_steplast, gamma_steplast,dc_steplast,vp_steplast,ds_steplast = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, finf,ainf,ginf,dcinf,vpinf,dsinf,-1,op_dir,date_time)
    beta_steplast, alpha_steplast, gamma_steplast, dc_steplast, vp_steplast, ds_steplast = PlasticityExperimentGauss(x_grid, [beta_array, alpha_arr, gamma_arr, dc_arr, vp_arr, ds_arr],
                              [r"$\beta$", r"$\alpha$", r"$\gamma$", r"$D_c$", r"$V_p$", r"$D_s$"], location,
                              l_scalesinf, [finf, ainf, ginf, dcinf, vpinf, dsinf], [1, 1, -1, 1, 1, 1], -1, op_dir, date_time)
    t_steplast = 29*60  # running for 10 secs
    b_factors.append(finf)
    a_factors.append(ainf)
    g_factors.append(ginf)
    time_steps.append(t_steplast)
    vp_factors.append(vpinf)
    dc_factors.append(dcinf)
    ds_factors.append(dsinf)
    # breakpoint()
    """
    chaning y_init to the last time_step value of step 2
    """
    new_y_init = data_mat[:, -1]
    model_params2 = [dc_steplast, ds_steplast, vp_steplast, Lamda_pc, Lamda_ps, alpha_steplast , beta_steplast, eta, omega,
                     gamma_steplast,
                     Jcin, Jsin, dx]

    t_evallast = np.arange(0, t_steplast, dt)
    t_rangelast = [0, t_steplast]
    soln_last = DynamicSimRun(model_params2, t_rangelast, t_evallast, new_y_init, max_step=100 * dt, method='RK45')
    print("Last finished at simulation time  = ", sim_time)

    """
    Concatinating the results of all steps
    """
    total_tps = np.concatenate((total_tps, (soln_last.t + sim_time)))
    sim_time += t_steplast
    print("total simulation time = ", sim_time)
    data_mat = np.concatenate((data_mat, soln_last.y), axis=1)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    savesimsettings(3, type,time_steps, lo,"{0}/protocol_{1}.json".format(op_dir, date_time),beta_factors=b_factors,
                    gamma_factors=g_factors, dc_factors=dc_factors,vp_factors=vp_factors,ds_factors=ds_factors,alpha_factors=a_factors)


beta_array = beta_orig * np.ones(P_c_init.shape)
gamma_arr = gamma_orig*np.ones(P_c_init.shape)
dc_arr = D_c_orig*np.ones(P_c_init.shape)
vp_arr = V_p_orig*np.ones(P_c_init.shape)
ds_arr = D_s_orig*np.ones(P_c_init.shape)
alpha_arr = alpha_orig*np.ones(P_c_init.shape)
eta_arr = eta_orig*np.ones(P_c_init.shape)
# breakpoint()
loc = [50]
location = [int(l) for l in loc]
# x_span = 3
# x_span_dx = int(x_span )
patterson2010GluA1Plasticity(ds_arr,dc_arr,vp_arr,alpha_arr,beta_array,eta_arr,gamma_arr, location)
