import AMPA_model
from AMPA_model import RunSimGluA2, RunModelWithFile,RunSimGluA1

import json
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.integrate import solve_ivp
import shutil


subunit= AMPA_model.subunit
baseline_param_file = "./ModelParamsTemporal{}.json".format(subunit)
dt = 0.02
P_s_init, P_c_init, P_spine_init, SP_model1 = RunModelWithFile(baseline_param_file)
y0 = np.vstack((np.vstack((P_c_init, P_s_init)), P_spine_init))
L = 500.0
dx = SP_model1.dx
D_s_orig = SP_model1.D_s  # in uM^2/s
D_c_orig = SP_model1.D_c  # in uM^2/s
V_p_orig = SP_model1.V_p  # in uM/s
half_life_surf = SP_model1.half_life_surf  # in days
Lamda_ps = np.log(2) / (half_life_surf * 24 * 60 * 60);
half_life_int = SP_model1.half_life_int  # in days
Lamda_pc = np.log(2) / (half_life_int * 24 * 60 * 60);
alpha_orig = SP_model1.alpha
beta_orig = SP_model1.beta;
Jcin = SP_model1.Jcin;
Jsin = SP_model1.Jsin;
omega = SP_model1.omega;  # same as eta_s_max in steady state model file
eta_orig = SP_model1.eta;  # same as eta_s_zero in steady state model file
gamma_orig = SP_model1.gamma
y_init_orig = y0.T.flatten()  # np.zeros(y0.T.flatten().shape) #np.zeros(y0.T.flatten().shape) #
print("baseline parameters ",dir(SP_model1))
# breakpoint()
def Endocytosis(ps, alpha):
    return alpha * ps;


def Exocytosis(pc, beta):
    return pc * beta;


def Degrdation(p, lamda_p):
    return p * lamda_p


def SpineExchange(ps, pspine, eta, omega, gamma):
    return eta * ps * (omega - pspine) - gamma * pspine;


def Spine_rk(p_s, p_spine, eta, omega, gamma):
    f_vec = lambda p_s, p_spine: eta * p_s * (omega - p_spine) - gamma * p_spine
    # print(f_vec)
    k1 = f_vec(p_s, p_spine)
    k2 = f_vec(p_s + k1 * dt / 2.0, p_spine + k1 * dt / 2.0)
    k3 = f_vec(p_s + k2 * dt / 2.0, p_spine + k2 * dt / 2.0)
    k4 = f_vec(p_s + dt * k3, p_spine + k3 * dt)

    return (1 / 6) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def AMPATimeDynamics(t, y, Dc, Ds, Vp, lamdac, lamdas, alpha, beta, eta, omega, gamma, Jcin, Jsin, dx):
    """
   Differential equations for the 1-D coupled diffusion-advection-reaction-trapping equations.

   The ODEs are derived using the method of lines.
   """
    # The vectors pc, ps and pspineare interleaved in y.  We define
    # views of pc,ps and pspine by slicing y.
    if t % 2 == 0:
        print("time = ", t, )
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

    # based on the updated and rechecked discretization calculations done 22 Aug 2024
    dpcdt[0] = Endocytosis(ps[0], alpha[0]) - Exocytosis(pc[0], beta[0]) - Degrdation(pc[0], lamdac) \
               + Jcin / dx + (Dc[0]/dx**2 - Vp[0]/dx )*pc[1] - (Dc[0]/dx**2) * pc[0]
    dpcdt[1:-1] = Endocytosis(ps[1:-1], alpha[1:-1]) - Exocytosis(pc[1:-1], beta[1:-1]) - Degrdation(pc[1:-1], lamdac) \
                  + Dc[1:-1] * np.diff(pc, 2) / dx**2 - Vp[1:-1] * np.diff(pc, 1)[1:] / dx
    dpcdt[-1] = Endocytosis(ps[-1], alpha[-1]) - Exocytosis(pc[-1], beta[-1]) - Degrdation(pc[-1], lamdac) \
                + (Dc[-1]/dx**2) * pc[-2] - (Dc[-1]/dx**2 + Vp[-1]**2/Dc[-1] + Vp[-1]/dx)*pc[-1]

    dpsdt[0] = Exocytosis(pc[0], beta[0]) - Endocytosis(ps[0], alpha[0]) - Degrdation(ps[0], lamdas) \
               - SpineExchange(ps[0], pspine[0], eta[0], omega,gamma[0]) + (Ds[0]/dx**2)*ps[1] - \
                (Ds[0]/dx**2)*ps[0] + Jsin/dx
    dpsdt[1:-1] = Exocytosis(pc[1:-1], beta[1:-1]) - Endocytosis(ps[1:-1], alpha[1:-1]) - Degrdation(ps[1:-1], lamdas) \
                  + Ds[1:-1] * np.diff(ps, 2) / dx**2 - SpineExchange(ps[1:-1], pspine[1:-1], eta[1:-1], omega, gamma[1:-1])
    dpsdt[-1] = Exocytosis(pc[-1], beta[-1]) - Endocytosis(ps[-1], alpha[-1]) - Degrdation(ps[-1], lamdas) \
                + (Ds[-1]/ dx**2) * ps[-2] - (Ds[-1]/dx**2)*ps[-1] - SpineExchange(ps[-1], pspine[-1], eta[-1], omega, gamma[-1])

    dpspinedt[0] = SpineExchange(ps[0], pspine[0], eta[0], omega, gamma[0])
    dpspinedt[1:-1] = SpineExchange(ps[1:-1], pspine[1:-1], eta[1:-1], omega, gamma[1:-1])  # Spine_rk(ps,pspine,eta,omega,gamma)
    dpspinedt[-1] = SpineExchange(ps[-1], pspine[-1], eta[-1], omega, gamma[-1])

    return dydt

# def Inerpparam(t,param):

def AdaptiveAMPATimeDynamics(t, y, Dc, Ds, Vp, lamdac, lamdas, alpha, beta_0,beta_mat, eta, omega, gamma, Jcin, Jsin, dx):
    # Getting the current value of exocytosis rate
    beta = beta_0*beta_mat[:,int(t/dt)]/beta_mat[:,0]

    # calculating the time derivative
    dydt = AMPATimeDynamics(t, y, Dc, Ds, Vp, lamdac, lamdas, alpha, beta, eta, omega, gamma, Jcin, Jsin, dx)
    return dydt
def DynamicSimRun(model_params, t_range, t_eval, y_init, method="LSODA", dense_op=True, vectorize=True, lband=2,
                  uband=2, rtol=1e-5, atol=1e-11, max_step=100):
    soln = solve_ivp(AMPATimeDynamics,
                     t_range,
                     y_init,
                     args=(model_params),
                     method=method,
                     dense_output=dense_op,
                     vectorize=vectorize,
                     lband=lband,
                     uband=uband,
                     rtol=rtol,
                     atol=atol,
                     max_step=max_step,
                     t_eval=t_eval)
    return soln

def AdaptiveDynamicsSimRun(model_params, t_range, t_eval, y_init, method="LSODA", dense_op=True, vectorize=True, lband=2,
                  uband=2, rtol=1e-5, atol=1e-11, max_step=100):
    # Passes the adptive version of temporal integrator
    soln = solve_ivp(AdaptiveAMPATimeDynamics, t_range, y_init, args=(model_params), \
                     method=method, dense_output=dense_op, vectorize=vectorize, lband=lband, uband=uband, rtol=rtol,
                     atol=atol, max_step=max_step, t_eval=t_eval)
    return soln
def SaveFigures(filename,ext_list = [".png",".svg",".pdf"],dpi=300):
    """
        function to save figures
        required arguments:
            filename
    """
    for ext in ext_list:
        plt.savefig(filename+ext,dpi=dpi)


def saveoutput(op_dir, date_time, soln, tps, percent, orig_params_file):
    p_c, p_s, p_sp = GetProteindata(soln)
    with open('{0}/PC_{1}_{2}_percent.npy'.format(op_dir, date_time, percent), 'wb') as f:
        np.save(f, p_c)
    f.close()
    with open('{0}/PS_{1}_{2}_percent.npy'.format(op_dir, date_time, percent), 'wb') as f:
        np.save(f, p_s)
    f.close()
    with open('{0}/PSPINE_{1}_{2}_percent.npy'.format(op_dir, date_time, percent), 'wb') as f:
        np.save(f, p_sp)
    f.close()
    t_points = tps
    with open('{0}/timepoints_{1}_{2}_percent.npy'.format(op_dir, date_time, percent), 'wb') as f:
        np.save(f, t_points)
    f.close()
    shutil.copyfile(orig_params_file, '{0}/baseline_parameters_{1}.json'.format(op_dir, date_time))

def BC(ya, yb):
    return np.array([ya[1] + Jsin / D_s_orig, yb[1], ya[3] - V_p_orig * ya[2] / D_c_orig + Jcin / D_c_orig, yb[3] - V_p_orig * yb[2] / D_c_orig])


def CheckBC(p_s, p_c):
    dps_a = np.diff(p_s, 1)[0]
    dpc_a = np.diff(p_c, 1)[0]
    pc_a = p_c[0]
    ps_a = p_s[0]
    ya = [ps_a, dps_a, pc_a, dpc_a]
    dps_b = np.diff(p_s, 1)[-1]
    dpc_b = np.diff(p_c, 1)[-1]
    pc_b = p_c[-1]
    ps_b = p_s[-1]
    yb = [ps_b, dps_b, pc_b, dpc_b]
    print("ya, yb = ", ya, yb)
    bcs = BC(ya, yb)
    return bcs

def GetProteindata(soln):
    pc_t = soln[::3]
    ps_t = soln[1::3]
    p_sp_t = soln[2::3]
    return pc_t, ps_t, p_sp_t

def savesimsettings(num_steps,type, time_steps, protocol_file,
                    dc_factors=[],ds_factors=[],vp_factors=[],
                    alpha_factors=[],beta_factors=[],
                    eta_factors=[],gamma_factors=[],
                    locations =[],unstim_locations = []):
    # assert num_steps == len(time_steps) - 1
    # assert num_steps == len(factors) - 1
    protocol_details = {}
    protocol_details["type"] = type
    protocol_details["num_steps"] = num_steps
    protocol_details["time_steps"] = time_steps
    protocol_details["stim locations"] = locations
    protocol_details["unstim locations"] = locations
    # protocol_details["x_span"] = x_span

    protocol_details["dc_factors"] = dc_factors
    protocol_details["ds_factors"] = ds_factors
    protocol_details["vp_factors"] = vp_factors
    protocol_details["alpha_factors"] = alpha_factors
    protocol_details["beta_factors"] = beta_factors
    protocol_details["alpha_factors"] = alpha_factors
    protocol_details["eta_factors"] = eta_factors
    protocol_details["gamma_factors"] = gamma_factors

    with open (protocol_file, 'a') as fp:
        json.dump(protocol_details,fp)
    fp.close()


def GlobalPlasticityExperiment(param,uod_factor):
    param_updated = param.copy()
    param_updated *= uod_factor
    return param_updated

def PlasticityExperimentGauss(x_grid,params,param_names, locns, sigmas, factors,up_or_down,step_num,op_dir,dt):
    new_pa = params.copy()
    fig,ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    plt.yscale("log")
    for pdx,p in enumerate(params):
        if sigmas[pdx] == dx:
            # print("modifying ",param_names[pdx],((factors[pdx])**up_or_down[pdx]))
            # for l1 in locns:
            new_pa[pdx] = Singledxchange(params[pdx],locns,factors[pdx],up_or_down[pdx])
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
    SaveFigures("{0}/fig_protocol_{1}_step_{2}".format(op_dir, dt, step_num), dpi=300)
    plt.show()
    return new_pa

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

def Singledxchange(pa,locns,factor,uod):
    p = pa.copy()
    for l1 in locns:
        p[int(l1/dx)]  *= ((factor)**uod)
    return p