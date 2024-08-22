from datetime import datetime
from GetPublishedData import oh_data_stim_sa,oh_data_unstim_sa
from TemporalIntegration import *

def Oh2051MultiSpinePlasticity(ds,dc,vp,alpha,beta,eta,gamma,stim_locs,unstim_locs):
    """
    step zero, we simulate the model in steady state for 30 secs
    step First, we step-increase the exocytosis rate (beta) to f1 folds and spine uptake rate (eta) to e1 folds in 10 mins
    step last, we change back the parameters to basal level and run simulation for 50 mins
    we integrate for a total of 30 mins to see the GluA1 dynamics
    """
    type = "Mutti_spine_Stim"
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/Time-dependent/" + date_time
    os.makedirs(op_dir, exist_ok=True)

    sim_time = 0
    time_steps = [0]
    b_factors = []
    a_factors = []
    e_factors = []
    num_steps =0

    """
    Step 0
    """
    b0 = 1  # increase by a factor of 1
    e0 = 1
    a0 = 1
    t_step0 = 30  # running for t_step0 secs
    b_factors.append(b0)
    a_factors.append(a0)
    # g_factors.append(g0)
    e_factors.append(e0)
    time_steps.append(t_step0)
    beta_step0 = beta*b0
    eta_step0 = eta*e0
    alpha_step0 = alpha*a0
    model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha_step0, beta_step0, eta_step0, omega, gamma,
                    Jcin, Jsin, dx]
    num_steps += 1
    # plt.plot(x_grid, beta_step0 / beta_step0[0])
    # plt.show()
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln = DynamicSimRun(model_params, t_range, t_eval, y_init_orig, max_step=100 * dt, method='RK45')
    data_mat = soln.y
    total_tps = soln.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)
    # beta_profile = np.ones(soln0.t.shape) * beta_step0[location]
    # eta_profile = np.ones(soln0.t.shape) * eta_step0[location]

    """
    Step 1
    """
    t_current = 0
    eta_stim,eta_unstim = oh_data_stim_sa(),oh_data_unstim_sa()
    for s,u in zip(eta_stim,eta_unstim):
        sim_time = s[0]
        eta_up = 1+ (s[1]-100)/100
        eta_down = 1+(u[1] - 100)/100
        print(eta_up,eta_down)
        l_scale = [dx]
        eta_step = eta
        e_factors.append([eta_up,eta_down])
        eta_step, = PlasticityExperimentGauss(SP_model1.x_grid, [eta_step],[r"$\eta$"], stim_locs,l_scale, [eta_up],
                                              [1], 2*num_steps, op_dir, date_time)
        eta_step, = PlasticityExperimentGauss(SP_model1.x_grid, [eta_step], [r"$\eta$"], unstim_locs, l_scale, [eta_down],
                                              [1], 2*num_steps+1, op_dir, date_time)
        t_end = (sim_time-t_current)*60
        t_range = [0,t_end]

        model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha, beta, eta_step, omega, gamma,
                        Jcin, Jsin, dx]

        t_eval = np.arange(0, t_end, dt)
        new_y_init = data_mat[:, -1]
        soln = DynamicSimRun(model_params, t_range, t_eval, new_y_init, max_step=100 * dt, method='RK45')
        total_tps = np.concatenate((total_tps, (soln.t + (t_current)*60)))
        # sim_time += t_step1
        data_mat = np.concatenate((data_mat, soln.y), axis=1)
        t_current = sim_time
        num_steps += 1
        print("Current simulation time  = ", t_current*60)
    print("total simulation time = ", t_current)
    total_tps = np.concatenate((total_tps, (soln.t + (t_current)*60)))
    data_mat = np.concatenate((data_mat, soln.y), axis=1)
    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    savesimsettings(num_steps, type, time_steps, "{0}/protocol_{1}.json".format(op_dir, date_time),
                                    eta_factors=e_factors,
                                    locations=stim_locs,
                                    unstim_locations=unstim_locs)
    # b1 = 10
    # e1 = 1.13
    # eu1 = 1 -.07
    # a1 = 5
    # l_scales1 = [3,dx,5]
    # new_y_init = data_mat[:, -1]  # Bleach(data_mat[:, -1],lo,x_sdx)
    # beta_step1 , eta_step1,alpha_step1 =  PlasticityExperimentGauss(SP_model1.x_grid, [beta, eta,alpha],[r"$\beta$", r"$\eta$",r"$\alpha$"], stim_locations,l_scales1, [b1,e1,a1], [1, 1,1], 1, op_dir, date_time)
    # eta_step1, = PlasticityExperimentGauss(SP_model1.x_grid, [eta_step1],[r"$\eta$"], unstim_location,l_scales1[1:2], [eu1], [1], 2, op_dir, date_time)
    # t_step1 = 5*60  # running for 50 secs
    # b_factors.append(b1)
    # a_factors.append(a1)
    # e_factors.append([e1,eu1])
    # time_steps.append(t_step1)
    # model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha_step1 , beta_step1, eta_step1, omega, gamma,
    #                 Jcin, Jsin, dx]
    # # breakpoint()
    # t_range = [0, t_step1]
    # t_eval = np.arange(0, t_step1, dt)
    # soln1 = DynamicSimRun(model_params, t_range, t_eval, new_y_init, max_step=100 * dt, method='RK45')
    # total_tps = np.concatenate((total_tps, (soln1.t + sim_time)))
    # sim_time += t_step1
    # data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    # print("Step 1 finished at simulation time  = ", sim_time)
    #
    # """
    # Step last
    # """
    # binf = 1
    # einf = 1.13
    # euinf = 1 -0.07
    # ainf = 1
    # beta_steplast = beta*binf
    # l_scaleslast = [dx]
    # eta_steplast, =  PlasticityExperimentGauss(SP_model1.x_grid, [eta],[ r"$\eta$"], stim_locations,l_scaleslast, [einf], [ 1], 3, op_dir, date_time)
    # eta_steplast, = PlasticityExperimentGauss(SP_model1.x_grid, [eta_steplast], [r"$\eta$"], unstim_location, l_scaleslast,
    #                                           [euinf], [1], 4, op_dir, date_time)
    # # breakpoint()
    # t_steplast = 35*60  # running for 10 secs
    # b_factors.append(binf)
    # e_factors.append(einf)
    # a_factors.append(ainf)
    # time_steps.append(t_steplast)
    #
    # """
    # chaning y_init to the last time_step value of step 2
    # """
    # new_y_init = data_mat[:, -1]
    # model_params = [dc, ds, vp, Lamda_pc, Lamda_ps, alpha, beta_steplast, eta_steplast, omega,
    #                  gamma,
    #                  Jcin, Jsin, dx]
    # # plt.plot(x_grid, beta_steplast / beta_steplast[0])
    # # plt.show()
    # t_evallast = np.arange(0, t_steplast, dt)
    # t_rangelast = [0, t_steplast]
    # soln_last = DynamicSimRun(model_params, t_rangelast, t_evallast, new_y_init, max_step=100 * dt, method='RK45')
    # print("Last finished at simulation time  = ", sim_time)
    #
    # """
    # Concatinating the results of all steps
    # """
    #
    # print("date and time:", date_time)
    # total_tps = np.concatenate((total_tps, (soln_last.t + sim_time)))
    # sim_time += t_steplast
    # print("total simulation time = ", sim_time)
    # data_mat = np.concatenate((data_mat, soln_last.y), axis=1)
    # saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)
    # savesimsettings(3, type,time_steps,  "{0}/protocol_{1}.json".format(op_dir, date_time),
    #                 beta_factors=b_factors,
    #                 eta_factors=e_factors,
    #                 alpha_factors=a_factors,
    #                 locations=stim_locations,
    #                 unstim_locations=unstim_location)
    #
    #
    # plt.close()
    print("Simulation complete and data saved!")


beta_array = beta_orig * np.ones(P_c_init.shape)
gamma_arr = gamma_orig*np.ones(P_c_init.shape)
dc_arr = D_c_orig*np.ones(P_c_init.shape)
vp_arr = V_p_orig*np.ones(P_c_init.shape)
ds_arr = D_s_orig*np.ones(P_c_init.shape)
alpha_arr = alpha_orig*np.ones(P_c_init.shape)
eta_arr = eta_orig*np.ones(P_c_init.shape)
locations = [250,251,252,253,255,256]
us_locations = [254]
Oh2051MultiSpinePlasticity(ds_arr,dc_arr,vp_arr,alpha_arr,beta_array,eta_arr,gamma_arr,locations,us_locations)
