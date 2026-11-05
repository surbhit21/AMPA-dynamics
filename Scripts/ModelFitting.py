#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:54:20 2022

@author: surbhitwagle
"""
from AMPA_model import *
import numpy as np
from lmfit import conf_interval, minimize,Minimizer, Parameters, Parameter, report_fit, printfuncs
import matplotlib.pyplot as plt
import numpy as np
import os
import PlotBinned 
from PlottingWidgetAMPA import *
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
from functools import reduce
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp, kruskal


class ModelFitting():
    def FitModel(x,data,rat,soma_rat,pars=[]):
        if pars == []:
            fit_paramas = Parameters()
            np.random.seed(2022)
            dc_min = 0.05
            dc_max = 10
            dc_init = np.random.uniform(dc_min,dc_max)
            ds_min = 0.1
            ds_max = 1.0
            ds_init = np.random.uniform(ds_min,ds_max)
            vp_min = 0.0
            vp_max =1.5
            vp_init = np.random.uniform(vp_min,vp_max)
            
            # breakpoint()
            print("dc_init = ",dc_init)
            print("ds_init = ",ds_init)
            print("vp_init = ",vp_init)
            #  parameter ranges
            
            #  parameters to fit
            fit_paramas.add('dc',dc_init,min=0)
            fit_paramas.add('ds',ds_init,min=0)
            fit_paramas.add('vp',vp_init,min=0)
            
            # rat = 1.1 # ratio between alpha and beta
            #fixed parameters
            
            fit_paramas.add('dx',scale,vary=False)
            fit_paramas.add('half_life_surf',float('inf'),vary=False)
            fit_paramas.add('half_life_int',1.95,vary=False)
            fit_paramas.add('alpha',1,vary=False)
            fit_paramas.add('beta',1/rat,vary=False)
            fit_paramas.add('Jsin',0.021/soma_rat,vary=False)
            fit_paramas.add('Jcin',0.021,vary=False)
            fit_paramas.add('eta_s_max',60,vary=False)
            fit_paramas.add('eta_s_zero',1/(15*60),vary=False)
            fit_paramas.add('gamma',1/(15*60),vary=False)
        else:
            fit_paramas = pars
        
        # breakpoint()
        mini = Minimizer(resudual,fit_paramas,fcn_kws={'x':x, 'data':data})
        
        # out1 = minimize(resudual,fit_paramas,method='Nelder',tol=1e-10,args=(x, data))
        # report_fit(out1.params)
        out2 = minimize(resudual,params=fit_paramas,method='leastsq',args=(x, data))
        report_fit(out2.params)
        breakpoint()
        # breakpoint()
        # ci, trace = conf_interval(mini, out2, sigmas=[1, 2], trace=True)
        # printfuncs.report_ci(ci)
        # breakpoint()
        
        return FittedCalculation(out2.params,x,data)
            