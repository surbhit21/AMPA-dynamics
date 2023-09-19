import numpy as np
import lmfit
from lmfit import  minimize, Parameters, report_fit
from scipy.optimize import curve_fit
from pathlib import Path

cmaps = {}

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))
COLORS = ["#005f73","#9b2226","#CA6702","#337357"]
COLORS_dict = {"spine":"#006b05","shaft":'#CA6702',"spine_s":"#006b05","spine_i":'#90ec7c',"shaft_s":"#005f73","shaft_i":'#CA6702'}

def plot_color_gradients(category, cmap_list):
    # Create figure and adjust figure height to number of colormaps
    nrows = len(cmap_list)
    figh = 0.35 + 0.15 + (nrows + (nrows - 1) * 0.1) * 0.22
    fig, axs = plt.subplots(nrows=nrows + 1, figsize=(6.4, figh))
    fig.subplots_adjust(top=1 - 0.35 / figh, bottom=0.15 / figh,
                        left=0.2, right=0.99)
    axs[0].set_title(f'{category} colormaps', fontsize=14)

    for ax, name in zip(axs, cmap_list):
        ax.imshow(gradient, aspect='auto', cmap=mpl.colormaps[name])
        ax.text(-0.01, 0.5, name, va='center', ha='right', fontsize=10,
                transform=ax.transAxes)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axs:
        ax.set_axis_off()

    # Save colormap list for later.
    cmaps[category] = cmap_list

def Area(radius):
    return np.pi*(radius**2)

def CreateFolderRecursive(folder):
    Path(folder).mkdir(parents=True, exist_ok=True)


def GetMatSum(mat,ax=0):
    return np.sum(mat,axis=ax)

def SortPunctas(ps,column=0):
    return ps[ps[:,column].argsort()]


def ChiSq(yd,y_fit,sigmas):
    nzs = np.nonzero(sigmas)
    # print(nzs)
    r_yd = np.take(yd,nzs)
    r_yf = np.take(y_fit,nzs)
    r_sgs = np.take(sigmas,nzs)
    # print("r_sgs = ",r_sgs)
    residuals = r_yd - r_yf
    chi_squ = np.sum((residuals/ r_sgs)**2)
    # print(residuals,r_sgs)
    return chi_squ


def BinnedSum( arr, bins, num_col=-1, name=None):
    # print(name)
    if len(arr.shape) == 2:
        rr, cc = arr.shape
        binned_sum = np.zeros((len(bins), cc))
        digitized = bins.searchsorted(arr[:, num_col])
        #             print(digitized,bins,arr[:,0])
        # breakpoint()
        digitized[0] = digitized[1]
        for c in range(0, cc):
            binned_sum[:, c] = np.bincount(digitized, weights=arr[:, c], minlength=len(bins))
        binned_sum[:, num_col] = bins
        return binned_sum[1:]
    else:
        print("quite not the shape", arr.shape)
        return np.zeros((len(bins), arr.shape[1]))

def BinnedAVG( arr, bins, num_col=-1, name=None):
    # print(name)
    if len(arr.shape) == 2:
        rr, cc = arr.shape
        binned_avg= np.zeros((len(bins), cc))
        digitized = bins.searchsorted(arr[:, num_col])
        #             print(digitized,bins,arr[:,0])

        digitized[0] = digitized[1]
        digi_counts = np.bincount(digitized)
        # breakpoint()
        for c1 in range(0, cc):
            binned_avg[:, c1] = (np.bincount(digitized, weights=arr[:, c1], minlength=len(bins))/np.array([max(1, v) for v in digi_counts]))[1:]
        binned_avg[:, num_col] = bins
        return binned_avg[1:]
    else:
        print("quite not the shape", arr.shape)
        return np.zeros((len(bins), arr.shape[1]))

def GetUniqueRows(mat):
    return np.unique(mat, axis=0)

def ExpFit2(xdata,ydata,sigmas,Fidx,Lidx,molecule):
    pars = Parameters()
    pars.add('amplitude',1,vary=False)
    pars.add('decay',1,min=0)
    mod = lmfit.models.ExponentialModel()
    out = mod.fit(ydata[Fidx:], pars, x=xdata[Fidx:])
    y_fit = exponential(xdata[Fidx:],-1.0/out.params['decay'])
    residuals = ydata[Fidx:]- y_fit[Fidx:]
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata[Fidx:]-np.mean(ydata[Fidx:]))**2)
    r_squared = 1 - (ss_res / ss_tot)
    chi_squ = np.sum((residuals/sigmas)**2)
    print("here",chi_squ)
    # breakpoint()
    return y_fit,r_squared,out.chisqr


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def normExponential(x, params):
    b = params['b'].value
    return np.exp(b * x)


def oneExponential(x, params):
    a = params['a'].value
    b = params['b'].value
    return a * np.exp(b * x)


def twoExponential(x, params):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    return a * np.exp(b * x) + c * np.exp(d * x)


def ExpFit(ftype, xdata, ydata, sigmas, Fidx, Lidx, molecule):
    """
        Fit a function to a given distribution

    """
    if ftype == "NormE":
        param_bounds = ([-np.inf], [np.inf])
        popt, pcov = curve_fit(normExponential, xdata[Fidx:], ydata[Fidx:], bounds=param_bounds, maxfev=5000)
        y_fit = normExponential(xdata, *popt)
    elif ftype == "1E":
        param_bounds = ([-np.inf, -np.inf], [+np.inf, 0])
        popt, pcov = curve_fit(oneExponential, xdata[Fidx:], ydata[Fidx:], bounds=param_bounds)
        y_fit = oneExponential(xdata, *popt)
    elif ftype == "2E":
        param_bounds = ([-np.inf, -np.inf, -np.inf, -np.inf], [+np.inf, 0, +np.inf, 0])
        popt, pcov = curve_fit(twoExponential, xdata[Fidx:], ydata[Fidx:], bounds=param_bounds)
        y_fit = twoExponential(xdata, *popt)
    else:
        raise NotImplementedError("ftype: {} not implemented, contact author or define it yourself".format(ftype))

    print("fitted " + ftype, popt)
    residuals = ydata[Fidx:] - y_fit[Fidx:]
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata[Fidx:] - np.mean(ydata[Fidx:])) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    chi_squ = ChiSq(ydata[Fidx:], y_fit[Fidx:], sigmas[Fidx:])
    print("chi-squared = ", chi_squ)
    return y_fit, r_squared, chi_squ


def ExpFitWithMinimize(ftype, xdata, ydata, sigmas, Fidx, Lidx, molecule):
    """
        Fit a function to a given distribution using lmfit minimize method

    """
    fit_paramas = Parameters()
    # np.random.seed(2022)
    exp_min = -2.0
    exp_max = 0
    pref_min = 0
    pref_max = 200

    if ftype == "NormE":
        b_init = np.random.uniform(exp_min, exp_max)
        fit_paramas.add('b', b_init, min=exp_min, max=exp_max)

        residuals = Residual(fit_paramas, normExponential, xdata[Fidx:], ydata[Fidx:])
        out2 = minimize(Residual, params=fit_paramas, method='leastsq',
                        args=(normExponential, xdata[Fidx:], ydata[Fidx:]))
        report_fit(out2.params)
        y_fit = normExponential(xdata[Fidx:], out2.params)

        # breakpoint()

    elif ftype == "1E":

        a_init = np.random.uniform(pref_min, pref_max)
        b_init = np.random.uniform(exp_min, exp_max)
        fit_paramas.add('a', a_init, min=pref_min, max=pref_max)
        fit_paramas.add('b', b_init, min=exp_min, max=exp_max)
        residuals = Residual(fit_paramas, oneExponential, xdata[Fidx:], ydata[Fidx:])
        out2 = minimize(Residual, params=fit_paramas, method='leastsq',
                        args=(oneExponential, xdata[Fidx:], ydata[Fidx:]))
        report_fit(out2.params)
        y_fit = oneExponential(xdata[Fidx:], out2.params)

    elif ftype == "2E":
        a_init = np.random.uniform(pref_min, pref_max)
        b_init = np.random.uniform(exp_min, exp_max)
        c_init = np.random.uniform(pref_min, pref_max)
        d_init = np.random.uniform(exp_min, exp_max)
        fit_paramas.add('a', a_init, min=pref_min, max=pref_max)
        fit_paramas.add('b', b_init, min=exp_min, max=exp_max)
        fit_paramas.add('c', c_init, min=pref_min, max=pref_max)
        fit_paramas.add('d', d_init, min=exp_min, max=exp_max)
        residuals = Residual(fit_paramas, twoExponential, xdata[Fidx:], ydata[Fidx:])
        out2 = minimize(Residual, params=fit_paramas, method='leastsq',
                        args=(twoExponential, xdata[Fidx:], ydata[Fidx:]))
        report_fit(out2.params)
        y_fit = twoExponential(xdata[Fidx:], out2.params)
    else:
        raise NotImplementedError("ftype: {} not implemented, contact author or define it yourself".format(ftype))

    # print("fitted "+ ftype, popt)

    chi_squ = out2.chisqr

    return y_fit, chi_squ


def Residual(paras, fun, x, data):
    expected_vals = fun(x, paras)
    res = expected_vals - data
    return res

def getminmax(arr,orig_min=np.Inf,orig_max=-np.Inf):
    if np.min(arr) < orig_min:
        orig_min = np.min(arr)
    if np.max(arr) > orig_max:
        orig_max = np.max(arr)
    return orig_min,orig_max


def exp_fit(x,a,b):
    return a*np.exp(-b*x)

def line_fit(x, m,c):
    return m*x + c
def R_seq(y_fit,y_orig):
    ss_res = ((y_orig-y_fit)**2).sum()
    ss_tot = ((y_orig-y_orig.mean())**2).sum()
    # print("in R_Seq =",ss_tot,ss_res)
    return 1 - (ss_res/ss_tot)