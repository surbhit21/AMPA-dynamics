#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
from datetime import datetime

# import matplotlib.pyplot as plt

# import numpy as np
from TemporalIntegration import *





def Tanaka2012GluA1Plasticity(ds,dc,vp,alpha,beta,eta,gamma):
    """
    step zero, we simulate the model in steady state for 30 secs
    step First, we step-increase the exocytosis rate (beta) to f1 folds and spine uptake rate (eta) to e1 folds in 10 mins
    step last, we change back the parameters to basal level and run simulation for 50 mins
    we integrate for a total of 30 mins to see the GluA1 dynamics
    """
    type = "G2_cLTP"
    sim_time = 0
    time_steps = [0]
    b_factors = [1]
    e_factors = [1]
    """
    Step 0
    """
    b0 = 1  # increase by a factor of 1
    e0 = 1
    t_step0 = 30  # running for t_step0 secs
    b_factors.append(b0)
    # g_factors.append(g0)
    e_factors.append(e0)
    time_steps.append(t_step0)
    beta_step0 = GlobalPlasticityExperiment(beta,b0)
    eta_step0 = GlobalPlasticityExperiment(eta, e0)

    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha, beta_step0, eta_step0, omega, gamma,
                    Jcin, Jsin, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, y_init_orig, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)

    """
    Step 1
    """
    b1 = 10
    e1 = 1.3
    new_y_init = data_mat[:, -1]  # Bleach(data_mat[:, -1],lo,x_sdx)
    beta_step1 = GlobalPlasticityExperiment(beta, b1)
    eta_step1 = GlobalPlasticityExperiment(eta, e1)
    t_step1 = 10*60  # running for 50 secs
    b_factors.append(b1)
    e_factors.append(e1)
    time_steps.append(t_step1)
    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha , beta_step1, eta_step1, omega, gamma,
                    Jcin, Jsin, dx]
    t_range = [0, t_step1]
    t_eval = np.arange(0, t_step1, dt)
    soln1 = DynamicSimRun(model_params, t_range, t_eval, new_y_init, max_step=100 * dt, method='RK45')
    total_tps = np.concatenate((total_tps, (soln1.t + sim_time)))
    sim_time += t_step1
    data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    print("Step 1 finished at simulation time  = ", sim_time)
    # breakpoint()

    """
    Step last
    """
    binf = 1
    einf = 1.3
    beta_steplast = GlobalPlasticityExperiment(beta, binf)
    eta_steplast = GlobalPlasticityExperiment(eta,einf)
    t_steplast = 50*60  # running for 10 secs
    b_factors.append(binf)
    e_factors.append(einf)
    time_steps.append(t_steplast)

    """
    chaning y_init to the last time_step value of step 2
    """
    new_y_init = data_mat[:, -1]
    model_params2 = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha, beta_steplast, eta_steplast, omega,
                     gamma,
                     Jcin, Jsin, dx]
    t_evallast = np.arange(0, t_steplast, dt)
    t_rangelast = [0, t_steplast]
    soln_last = DynamicSimRun(model_params2, t_rangelast, t_evallast, new_y_init, max_step=100 * dt, method='RK45')
    print("Last finished at simulation time  = ", sim_time)

    """
    Concatinating the results of all steps
    """
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/Time-dependent/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    total_tps = np.concatenate((total_tps, (soln_last.t + sim_time)))
    sim_time += t_steplast
    print("total simulation time = ", sim_time)
    data_mat = np.concatenate((data_mat, soln_last.y), axis=1)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    savesimsettings(1, type,time_steps,  "{0}/protocol_{1}.json".format(op_dir, date_time), beta_factors=b_factors,
                    eta_factors=e_factors)
    plt.close()
    print("Simulation complete and data saved!")


beta_array = beta_orig * np.ones(P_c_init.shape)
gamma_arr = gamma_orig*np.ones(P_c_init.shape)
dc_arr = D_c_orig*np.ones(P_c_init.shape)
vp_arr = V_p_orig*np.ones(P_c_init.shape)
ds_arr = D_s_orig*np.ones(P_c_init.shape)
alpha_arr = alpha_orig*np.ones(P_c_init.shape)
eta_arr = eta_orig*np.ones(P_c_init.shape)
# loc = 250
# location = int(loc / dx)
# x_span = 250
# x_span_dx = int(x_span / dx)
Tanaka2012GluA1Plasticity(ds_arr,dc_arr,vp_arr,alpha_arr,beta_array,eta_arr,gamma_arr)
