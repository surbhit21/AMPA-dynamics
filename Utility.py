import csv
# import itertools
# import json
# # from lmfit import conf_interval, minimize,Minimizer, Parameters, Parameter, report_fit, printfuncs
# # from matplotlib.cbook import boxplot_stats
import json
import numpy as np
import lmfit
from lmfit import  minimize, Parameters, report_fit
import os
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit
cmaps = {}

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))
COLORS = ["#005f73","#9b2226","#CA6702","#337357"]
COLORS_dict = {"spine":"#005f73",
               "shaft":'#43AA8B',
               "spine_s":"#F9C74F",
               "spine_i":'#CA6702',
               "shaft_s":"#F94144",
               "shaft_i":'#F8961E',
               "s_2_i_ratio":"#43AA8B",
               "spine_2_shaft_ratio":"#4D908E",
               "soma":"#277DA1"}


# spine to dendrite ratio calculation setting
conditions = [ "TTX", "Untreated"]#,"control" ]
root_inp = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/"
spine_cell_folder = {
    "GluA2":{
        "control":"{}GluA2/Control/stretches/retake".format(root_inp),
        "TTX":"{}GluA2/TTX/New_data/stretches/TTX".format(root_inp),
        "Untreated":"{}GluA2/TTX/New_data/stretches/PBS".format(root_inp)
    },
    "GluA1": {
        "control":"{}GluA1/new/Control/stretches".format(root_inp),
        # "TTX":"{}GluA2/TTX/New_data/stretches/TTX".format(root_inp),
        # "Untreated":"{}GluA2/TTX/New_data/stretches/PBS".format(root_inp)},
    }
}
spine_num_cells = {
    "GluA2":{
        "control":12,
        "TTX":13,
        "PBS":10
    },
    "GluA1": {
            "control":6,
            # "TTX":13,
            # "PBS":10
    }
}
spine_exclude_cells = {
    "GluA2":{
        "control":[2,3,4,10],
        "TTX":[4,5,6,8,10,11,12],
        "PBS":[1,2]
    },
    "GluA1": {
            "control":[5],
            # "TTX":13,
            # "PBS":10
    }
}
spine_channel_names = {
    "GluA2":["GFP", "surf-GluA2", "int-GluA2"],
    "GluA1": ["GFP", "surf-GluA1", "int-GluA1"]
}
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

def RedChisq(yd,y_fit,sigmas,dof):
    return ChiSq(yd,y_fit,sigmas)/dof
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


def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))
def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
    return _1gaussian(x,amp1,cen1,sigma1) + _1gaussian(x,amp2,cen2,sigma2)

class GluA2StoD():
    def __init__(self, DataFolder, projection_type,channel_names):
        self.df = DataFolder
        self.molecules = channel_names
        self.GFP_channel = 0
        self.surf_channel = 1
        self.int_channel = 2
        self.projection_type = projection_type
        self.f_names = ["all_data.csv"]

    def LoadData(self, num_cells, exclude_cells=[]):
        """
        Assumes a folder strucutre. Give the outermost folder which contains folder for each image
           df ==> cell_i ==> Rois ==> New ==> int-GLuA2/suef-GluA2/GFP
        """

        # fig,(ax0,ax1) = plt.subplots(1, 2,figsize=(20, 12), sharey=True)
        files = os.listdir(self.df)
        stat_names = ["area", "mean", "stddev", "min", "max", "intden", "median", "rawintden", "RRLength", "RRWidth"]
        conditions = ["sp_", "sp_bg_", "dend_", "dend_bg_"]
        int_df = []
        surf_df = []
        gfp_df = []
        int_data = {}
        surf_data = {}
        gfp_data = {}

        int_dend_data = {}
        int_spine_line_data = {}
        surf_dend_data = {}
        surf_spine_line_data = {}
        gfp_dend_data = {}
        gfp_spine_line_data = {}

        for fdx in range(1, num_cells + 1):
            # print(file)
            file = "cell_{}".format(fdx)
            if os.path.isdir(os.path.join(self.df, file)) and fdx not in exclude_cells:
                sub_folder = "/Rois/{}/".format(self.projection_type)
                measures = "/measures"
                data = "/data"
                file1 = file + sub_folder
                spine_file = "/Spine/Synapse_l.json"

                int_folder = os.path.join(self.df, file1 + self.molecules[self.int_channel] + measures);
                surf_folder = os.path.join(self.df, file1 + self.molecules[self.surf_channel] + measures)
                GFP_folder = os.path.join(self.df, file1 + self.molecules[self.GFP_channel] + measures)

                spine_json_file = os.path.join(self.df, file + spine_file)
                # print("spine_file = ",spine_json_file)
                dist_dict = self.GetSpineDistance(spine_json_file, fdx)
                k = 0
                for i in range(0, len(self.f_names)):
                    int_data, num_samples = self.ReadCSVFull(int_folder + "/" + self.f_names[i], int_data, stat_names,
                                                             fdx, distance_dict=dist_dict)
                    # int_data.append(self.ReadCSVFull(int_folder+"/"+self.f_names[i+1]).shape)
                    surf_data, num_samples = self.ReadCSVFull(surf_folder + "/" + self.f_names[i], surf_data,
                                                              stat_names, fdx, distance_dict=dist_dict)
                    # surf_data.append(self.ReadCSVFull(int_folder+"/"+self.f_names[i+1]).shape)
                    gfp_data, num_samples = self.ReadCSVFull(GFP_folder + "/" + self.f_names[i], gfp_data, stat_names,
                                                             fdx, distance_dict=dist_dict)

                    int_dend_data["dend_{}".format(fdx)] = self.ReadCSVFull(
                        os.path.join(self.df, file1 + self.molecules[self.int_channel] + data + "/Dendrite1.csv"),
                        p_or_m=0)
                    surf_dend_data["dend_{}".format(fdx)] = self.ReadCSVFull(
                        os.path.join(self.df, file1 + self.molecules[self.surf_channel] + data + "/Dendrite1.csv"),
                        p_or_m=0, distance_dict=dist_dict)
                    gfp_dend_data["dend_{}".format(fdx)] = self.ReadCSVFull(
                        os.path.join(self.df, file1 + self.molecules[self.GFP_channel] + data + "/Dendrite1.csv"),
                        p_or_m=0, distance_dict=dist_dict)
                    int_spine_line_data["cell_{}".format(fdx)] = {}
                    surf_spine_line_data["cell_{}".format(fdx)] = {}
                    gfp_spine_line_data["cell_{}".format(fdx)] = {}
                    for spx in range(1, num_samples + 1):
                        int_sp_i_line_data = self.ReadCSVFull(os.path.join(self.df, file1 + self.molecules[
                            self.int_channel] + data + "/spines/SP{}-line.csv".format(spx)), p_or_m=0,
                                                              distance_dict=dist_dict)
                        surf_sp_i_line_data = self.ReadCSVFull(os.path.join(self.df, file1 + self.molecules[
                            self.surf_channel] + data + "/spines/SP{}-line.csv".format(spx)), p_or_m=0,
                                                               distance_dict=dist_dict)
                        gfp_sp_i_line_data = self.ReadCSVFull(os.path.join(self.df, file1 + self.molecules[
                            self.GFP_channel] + data + "/spines/SP{}-line.csv".format(spx)), p_or_m=0,
                                                              distance_dict=dist_dict)

                        int_spine_line_data["cell_{}".format(fdx)]["sp_{}".format(spx)] = np.asarray(int_sp_i_line_data)
                        surf_spine_line_data["cell_{}".format(fdx)]["sp_{}".format(spx)] = np.asarray(
                            surf_sp_i_line_data)
                        gfp_spine_line_data["cell_{}".format(fdx)]["sp_{}".format(spx)] = np.asarray(gfp_sp_i_line_data)
                # breakpoint()
                # int_spine_line_data["cell_{}".format(fdx)] = np.asarray(int_spine_line_data["cell_{}".format(fdx)],dtype=object)
                # surf_spine_line_data["cell_{}".format(fdx)] = np.asarray(surf_spine_line_data["cell_{}".format(fdx)],dtype=object)
                # gfp_spine_line_data["cell_{}".format(fdx)] = np.asarray(gfp_spine_line_data["cell_{}".format(fdx)],dtype=object)
        col_names = [item_two + item_one for item_one in stat_names for item_two in conditions]
        col_names.append("DFO")
        int_df = pd.DataFrame.from_dict(int_data,
                                        orient='index', columns=col_names)
        surf_df = pd.DataFrame.from_dict(surf_data,
                                         orient='index', columns=col_names)
        gfp_df = pd.DataFrame.from_dict(gfp_data,
                                        orient='index', columns=col_names)

        return int_df, surf_df, gfp_df, int_dend_data, surf_dend_data, gfp_dend_data, int_spine_line_data, surf_spine_line_data, gfp_spine_line_data


    def SpineWiseDict(self, num_samples, data, data_dict, cols, cell_id, distance_dict):
        # rr,cc = data.shape
        # num_samples = rr//4
        # print(data)
        # data_dict = {}
        # print(distance_dict)
        # breakpoint()
        for i in range(1, num_samples + 1):
            data_dict["cell_{}_sp_{}".format(cell_id, i)] = []
            for jdx, col in enumerate(cols):
                # print("cell_{}_sp_{}".format(cell_id,i))
                # print(jdx,col,data['sp{}'.format(i)][jdx],data['sp{}-bg'.format(i)][jdx],data['dend{}'.format(i)][jdx],data['dend{}-bg'.format(i)][jdx])
                data_dict["cell_{}_sp_{}".format(cell_id, i)].append(float(data['sp{}'.format(i)][jdx]))
                data_dict["cell_{}_sp_{}".format(cell_id, i)].append(float(data['sp{}-bg'.format(i)][jdx]))
                data_dict["cell_{}_sp_{}".format(cell_id, i)].append(float(data['dend{}'.format(i)][jdx]))
                data_dict["cell_{}_sp_{}".format(cell_id, i)].append(float(data['dend{}-bg'.format(i)][jdx]))
            data_dict["cell_{}_sp_{}".format(cell_id, i)].append(float(distance_dict['sp{}'.format(i)]))
        return data_dict, num_samples

    def GetSpineDistance(self, fname, cell_id=None):
        distance_dict = {}
        with open(fname, "r") as file:
            # Load the contents of the file as JSON data
            fdata = json.load(file)
        for spdx, spdata in enumerate(fdata):
            # print("cell_{}_sp_{} distance = {}".format(cell_id,spdx+1,spdata["distance"]))
            distance_dict["sp{}".format(spdx + 1)] = spdata["distance"]
        return distance_dict

    def ReadCSVFull(self, filename, data_dict=None, cols=None, cell_id=None, p_or_m=1, distance_dict=None):
        """
       function to read CSV data files:
           argumentrs : filename : path of the .csv file
               data_dict: data dict to store the read data, used only in case p_or_m == 1
               cols: number of columns in the file, used only in case p_or_m == 1
               cell_id: id of the cell
               p_or_m: = 1 if it is a imageJ measure file =0 if it is a imagJ plot profile file

       """

        with open(filename) as csvDataFile:
            # read file as csv file
            csvReader = csv.reader(csvDataFile)
            # loop over rows
            num_row = -1
            if p_or_m == 1:
                csv_data = {}
                for row in csvReader:
                    # add cell [0] to list of dates
                    # print(row[1].split(':')[-1])
                    num_row += 1
                    key = row[1].split(':')[-1]
                    csv_data[key] = row[-1 * len(cols):]
                # print(csv_data)
                return self.SpineWiseDict(num_row // 4, csv_data, data_dict, cols, cell_id, distance_dict)
            elif p_or_m == 0:
                csv_data = []
                for row in csvReader:
                    csv_data.append(row)
                data = np.asarray(csv_data[1:]).astype("float")
                # print(data[0,-2:])
                # print(data[:,1:])
                return data[:, 1:]
            else:
                raise ("some other type file sent, not proccessible")
