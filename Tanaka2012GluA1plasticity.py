#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
from datetime import datetime

import matplotlib.pyplot as plt

# import numpy as np
from TemporalIntegration import *


def PlasticityExperiment(beta,gamma, lo, x_sdx, up_or_down_factor,g_factor):
    beta_updated = beta.copy()
    gamaa_updated = gamma.copy()
    print(up_or_down_factor,g_factor)
    beta_updated *= up_or_down_factor
    gamaa_updated*= g_factor
    return beta_updated, gamaa_updated


# def Bleach(y, lo, x_sdx):
#     pc_l, ps_l, psp_l = GetProteindata(y)
#     ps_l[lo - x_span_dx:lo + x_span_dx] = 0
#     psp_l[lo - x_span_dx:lo + x_span_dx] = 0
#     return np.vstack((np.vstack((pc_l, ps_l)), psp_l)).T.flatten()


def Tanaka2012GluA1Plasticity(ds,dc,vp,alpha,beta,eta,gamma, lo, x_sdx):
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
    g_factors = [1]
    """
    Step 0
    """
    f0 = 1  # increase by a factor of 1
    g0 = 1
    t_step0 = 30  # running for t_step0 secs
    factors.append(f0)
    g_factors.append(g0)
    time_steps.append(t_step0)
    beta_step0,gamma_step0 = PlasticityExperiment(beta,gamma, lo, x_sdx, f0, g0)

    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha, beta_step0, eta, omega, gamma_step0,
                    Jcin, Jsin, dx]
    # plt.plot(x_grid, beta_step0 / beta_step0[0])
    # plt.show()
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, y_init_orig, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)
    beta_profile = np.ones(soln0.t.shape) * beta_step0[location]
    gamma_profile = np.ones(soln0.t.shape) * gamma_step0[location]
    """
    Step 1
    """
    f1 = 10
    g1 = 1/1.3
    new_y_init = data_mat[:, -1]  # Bleach(data_mat[:, -1],lo,x_sdx)
    beta_step1,gamma_step1 = PlasticityExperiment(beta,gamma, lo, x_sdx, f1,g1)
    t_step1 = 10*60  # running for 50 secs
    factors.append(f1)
    g_factors.append(g1)
    time_steps.append(t_step1)
    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha , beta_step1, eta, omega, gamma_step1,
                    Jcin, Jsin, dx]
    # plt.plot(x_grid, beta_step1 / beta_step1[0])
    # plt.show()
    t_range = [0, t_step1]
    t_eval = np.arange(0, t_step1, dt)
    soln1 = DynamicSimRun(model_params, t_range, t_eval, new_y_init, max_step=100 * dt, method='RK45')
    total_tps = np.concatenate((total_tps, (soln1.t + sim_time)))
    sim_time += t_step1
    data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    print("Step 1 finished at simulation time  = ", sim_time)
    # saveoutput(op_dir, date_time, data_mat, total_tps, 1, "./ModelParams.json")
    beta_profile = np.concatenate((beta_profile, np.ones(soln1.t.shape) * beta_step1[location]))
    gamma_profile = np.concatenate((gamma_profile, np.ones(soln1.t.shape) * gamma_step1[location]))

    """
    Step last
    """
    finf = 1
    ginf = 1/1.3
    beta_steplast, gamma_steplast = PlasticityExperiment(beta,gamma, lo, x_sdx, finf,ginf)
    t_steplast = 20*60  # running for 10 secs
    factors.append(finf)
    g_factors.append(ginf)
    time_steps.append(t_steplast)
    # breakpoint()
    """
    chaning y_init to the last time_step value of step 2
    """
    new_y_init = data_mat[:, -1]
    model_params2 = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha, beta_steplast, eta, omega,
                     gamma_steplast,
                     Jcin, Jsin, dx]
    # plt.plot(x_grid, beta_steplast / beta_steplast[0])
    # plt.show()
    t_evallast = np.arange(0, t_steplast, dt)
    t_rangelast = [0, t_steplast]
    soln_last = DynamicSimRun(model_params2, t_rangelast, t_evallast, new_y_init, max_step=100 * dt, method='RK45')
    print("Last finished at simulation time  = ", sim_time)
    # saveoutput(op_dir, date_time, data_mat, total_tps, n, "./ModelParams.json")
    beta_profile = np.concatenate((beta_profile, np.ones(soln_last.t.shape) * beta_steplast[location]))
    gamma_profile = np.concatenate((gamma_profile, np.ones(soln_last.t.shape) * gamma_steplast[location]))

    """
    Concatinating the results of all steps
    """
    total_tps = np.concatenate((total_tps, (soln_last.t + sim_time)))
    sim_time += t_steplast
    print("total simulation time = ", sim_time)
    data_mat = np.concatenate((data_mat, soln_last.y), axis=1)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, "./ModelParamsTemporal.json")
    breakpoint()
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    ax.plot((total_tps-t_step0)/60, beta_profile/beta_profile[0],label=r"$\beta$")
    ax.plot((total_tps-t_step0)/60, gamma_profile/gamma_profile[0],label=r"$\gamma$")
    ax.spines[["top","right"]].set_visible(False)
    # plt.xlabel("Simulation time in mins")
    # plt.ylabel("Fold change")
    plt.legend(frameon=False)
    # plt.legend()
    # breakpoint()
    plt.tight_layout()
    plt.savefig("{0}/fig_protocol_{1}.pdf".format(op_dir, date_time),dpi=300)
    plt.savefig("{0}/fig_protocol_{1}.png".format(op_dir, date_time), dpi=300)
    plt.savefig("{0}/fig_protocol_{1}.eps".format(op_dir, date_time), dpi=300)
    savesimsettings(3, time_steps, [loc], x_span, "{0}/protocol_{1}.json".format(op_dir, date_time),beta_factors=factors,gamma_factors=g_factors)
    plt.show()


beta_array = beta_orig * np.ones(P_c_init.shape)
gamma_arr = gamma_orig*np.ones(P_c_init.shape)
dc_arr = D_c_orig*np.ones(P_c_init.shape)
vp_arr = V_p_orig*np.ones(P_c_init.shape)
ds_arr = D_s_orig*np.ones(P_c_init.shape)
alpha_arr = alpha_orig*np.ones(P_c_init.shape)
eta_arr = eta_orig*np.ones(P_c_init.shape)
loc = 250
location = int(loc / dx)
x_span = 250
x_span_dx = int(x_span / dx)
Tanaka2012GluA1Plasticity(ds_arr,dc_arr,vp_arr,alpha_arr,beta_array,eta_arr,gamma_arr, location, x_span_dx)
