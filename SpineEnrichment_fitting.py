import sys

# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1,
                '/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Mo/Python-code/ampa-dynamics')

# import file
from AMPA_model import *
import csv
import json
from lmfit import conf_interval, minimize, Minimizer, Parameters, Parameter, report_fit, printfuncs
import matplotlib.pyplot as plt
import numpy as np
import os
import PlotBinned
from PlottingWidgetAMPA import *
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
from functools import reduce
from scipy.optimize import curve_fit
# from scipy.stats import ks_2samp, kruskal
#
# # import scikit_posthocs as sps
# import seaborn as sns
# import SNSPlottingWidget as SNSP
from Utility import BinnedAVG

scale_GluA2_protein = 0.24034
Lengths = np.array([10,20,30,40,50,60])
bin_size = 5
bins = np.arange(0, Lengths.max(), bin_size)


def ExpFit(xdata, ydata):
    param_bounds = ([0, 0], [np.inf, np.inf])
    popt, pcov = curve_fit(GetAMPADistribution, xdata, ydata, bounds=param_bounds)
    print(popt)
    y_fit = exponential(xdata, *popt)
    residuals = ydata - y_fit
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    return y_fit, r_squared


def FitModel(x, data, sigmas, rat,soma_rat,length_to_fit, pars=[]):
    if pars == []:
        fit_paramas = Parameters()
        # np.random.seed(2022)

        #  parameter ranges
        eta_min = -6.
        eta_max = -3.
        eta_init = np.random.uniform(eta_min, eta_max)
        #  parameters to fit
        fit_paramas.add('dc', 0.7, vary=False)
        fit_paramas.add('ds', 0.7, vary=False)
        fit_paramas.add('vp', 1e-5,vary=False)

        # rat = 1.1 # ratio between alpha and beta
        # fixed parameters

        fit_paramas.add('dx', 0.24, vary=False)
        fit_paramas.add('half_life_surf', float('inf'), vary=False)
        fit_paramas.add('half_life_int', 3.12, vary=False)
        fit_paramas.add('alpha', 1.8e-4, vary=False)
        fit_paramas.add('beta', 1.8e-4 * rat, vary=False)
        fit_paramas.add('Jsin', 0.21 / soma_rat, vary=False)
        fit_paramas.add('Jcin', 0.21, vary=False)
        fit_paramas.add('omega', 60, vary=False)
        fit_paramas.add('eta', eta_init, min=eta_min,max=eta_max)
        fit_paramas.add('gamma', 1 / (43), vary=False)
        fit_paramas.add('length_to_fit', length_to_fit, vary=False)
    else:
        fit_paramas = pars

    # breakpoint()
    mini = Minimizer(resudual, fit_paramas, fcn_kws={'x': x, 'data': data})

    # out1 = minimize(resudual,fit_paramas,method='Nelder',tol=1e-10,args=(x, data))
    # report_fit(out1.params)
    out2 = mini.minimize()
    report_fit(out2)
    # breakpoint()
    # breakpoint()
    # ci, trace = conf_interval(mini, out2, sigmas=[1, 2], trace=True)
    # printfuncs.report_ci(ci)
    # breakpoint()

    return FittedCalculation(out2.params, x, data, sigmas, mini, out2)


def resudual(paras, x=None, data=None, ):
    # pc_needed = pc_needed/(pc_binned.sum()*delta_x)
    # breakpoint()
    SE_needed = GetRequiredDist(paras, x, data)

    resd = SE_needed - data[0]
    return resd  # resd.flatten()


def GetSlidingWindowMean(data, window_len, mode='same'):
    try:
        conv_window = np.ones(window_len) * (1 / window_len)
        sw_mean = np.convolve(data, conv_window, mode=mode)
        return sw_mean
    except:
        print("exception")


def GetSlidingWindowMeanMatrix(data, window_len, mode='same'):
    if len(data.shape) != 2:
        return ("data is not matrix ")
    print("here")
    op_matrix = []
    # op_matrix = np.ones((data.shape[0],op.shape[0]))

    for d in data:
        # breakpoint()
        op_matrix.append(GetSlidingWindowMean(d, window_len, mode))
    op_matrix = np.asarray(op_matrix)
    return op_matrix


def GetRequiredDist(paras, x, data):
    x1, SE= GetParamAndModelDist(paras)
    length_to_fit = paras['length_to_fit'].value
    ltf_dx = int(length_to_fit / paras['dx'].value)
    x1 = x1[0:ltf_dx]
    SE_model = SE[0:ltf_dx]
    # G2DA1 = GluA2DataAnalysis("/")
    # print(data.shape)
    # binning the model distribution in the same size
    # breakpoint()
    SE_binned = BinnedAVG(np.column_stack((x1, SE_model)), bins, 0)[1:data.shape[0] + 1, 1]
    # pc_binned = BinnedSum(np.column_stack((x1, pc_model)), bins, 0)[1:data.shape[1] + 1, 1]
    # taking the first N bins
    # ps_needed = ps_binned[0:x.shape[0]]
    # # ps_needed = ps_needed/(ps_binned.sum()*delta_x)
    # pc_needed = pc_binned[0:x.shape[0]]

    # normalizing with the first bin / same as soma
    # ps_needed = ps_needed / ps_needed[0]
    # pc_needed = pc_needed / pc_needed[0]
    return SE_binned


def GetParamAndModelDist(paras):
    # reading parameters to fit
    D_c = paras['dc'].value
    D_s = paras['ds'].value

    #  reading fixed parameters
    V_p = paras['vp'].value
    # print("*",end="")
    # print("Dc = ",D_c,"Ds = ",D_s,"Vp = ",V_p)
    delta_x = paras['dx'].value
    half_life_surf = paras['half_life_surf'].value
    half_life_int = paras['half_life_int'].value
    alpha = paras['alpha'].value
    beta = paras['beta'].value
    Jsin = paras['Jsin'].value
    Jcin = paras['Jcin'].value
    omega = paras['omega'].value
    eta = 10**paras['eta'].value
    gamma = paras['gamma'].value
    # breakpoint()
    # return model distribution
    return RunModel(D_s, D_c, V_p, half_life_surf, half_life_int, alpha, beta, Jsin, Jcin, omega, eta, gamma,
                    delta_x)


def RunModel(D_s, D_c, V_p, half_life_surf, half_life_int, alpha, beta, Jsin, Jcin, omega, eta, gamma,
             delta_x):
    SP_model1 = DendriteWithStochasticSpinesConstantV(D_s, D_c, V_p, half_life_surf, half_life_int, alpha, beta, Jsin,
                                                      Jcin, omega, eta, gamma, delta_x);
    ps_dist, pc_dist = SP_model1.solveNumerical()
    ps_spine = SP_model1.omega * (1 / (1 + (SP_model1.gamma / (SP_model1.eta * ps_dist))))
    x1 = SP_model1.x_grid
    SE = ps_spine/ps_dist
    print(pc_dist[1040:1080],ps_spine[1040:1080],ps_dist[1040:1080],SE[1040:1080],SP_model1.eta)
    # ps_sum, pc_sum, ps_ana, pc_ana = SP_model1.IntegralBC(ps_dist, pc_dist)
    # returning sum normalized distribution
    return x1,SE


def FittedCalculation(paras, x, data, sigmas, mini, out2):
    x1, SE_dist = GetParamAndModelDist(paras)
    SE_needed = GetRequiredDist(paras, x, data)
    # GetParamAndModelDist
    delta_x = paras['dx'].value
    x_n = int(np.ceil(x[-1] / delta_x))
    SE_chi2 = ChiSq(data[0], SE_needed, sigmas[0])
    return x1[0:x_n],SE_dist[0:x_n],SE_chi2,paras,mini,out2


def R_seq(ydata, y_fit):
    residuals = ydata - y_fit
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def ChiSq(yd, y_fit, sigmas):
    nzs = np.nonzero(sigmas)
    # print(nzs)
    r_yd = np.take(yd, nzs)
    r_yf = np.take(y_fit, nzs)
    r_sgs = np.take(sigmas, nzs)
    residuals = r_yd - r_yf
    chi_squ = np.sum((residuals / r_sgs) ** 2)
    print(residuals, r_sgs)
    return chi_squ


def plotFittedModelWithData(x1, mean_surf, mean_int, fit_surf, fit_int):
    fig, ax = plt.subplot(nrows=2, ncols=1)
    ax[0].plot

