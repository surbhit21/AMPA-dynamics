import csv
import seaborn as sns
import matplotlib.pyplot as plt
import numpy
import matplotlib
matplotlib.use("Qt5Agg")
import numpy as np
import os
import pandas as pd
from  SNSPlottingWidget import SNSPlottingWidget
from Utility import COLORS_dict
import scipy as sp
root_folder = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/"
Patterson_2010_4c_file = os.path.join(root_folder,"Patterson_2010_4C.csv")
Patterson_2010_1g_file = os.path.join(root_folder,"Patterson_2010_1G.csv")
Patterson_2010_1b_file = os.path.join(root_folder,"Patterson_2010_1B.csv")
Tanaka_2012_S3a_PSD_file = os.path.join(root_folder,"Tanaka_2012_S3A_PSLM.csv")
Tanaka_2012_S3a_Dend_file = os.path.join(root_folder,"Tanaka_2012_S3A_NON_PSLM.csv")
Tanaka_2012_exo_data = os.path.join(root_folder,"Tanaka_2012_exo_NPLSM.csv")

Graves_2021_glua1 = os.path.join(root_folder,"Graves_2021_5DGreen.csv")
Graves_2021_dsred = os.path.join(root_folder,"Graves_2021_5DRed.csv")

my_pal1= {"Control": COLORS_dict["spine"],"cLTP": COLORS_dict["shaft"]}
Clavet_Fournier_2023_4_folder = os.path.join(root_folder,"K_Willig/Fig4/")
Clavet_Fournier_2023_4_fnames = [
    'Fig4_Control_0min',
    'Fig4_Control+ 30 min',
    'Fig4_Control+ 1h',
    'Fig4_Control+ 2h',
    'Fig4_cLTP_0 min',
    'Fig4_cLTP+ 30 min',
    'Fig4_cLTP+ 60 min',
    'Fig4_cLTP + 2h'
]
def ReadCSVFull(filename,head_rows = 0):
    """
   function to read CSV data files:
       arguments : filename : path of the .csv file
   """
    data_mat = []
    with open(filename) as csvDataFile:
        # read file as csv file
        csvReader = csv.reader(csvDataFile)
        # loop over rows
        num_rows = 0
        for row in csvReader:
            # add cell [0] to list of dates
            if row != []:
                num_rows += 1
                row_data = [get_valid_float(x) for x in row]
                data_mat.append(row_data)
        assert num_rows%2 == 0
        data_mat = np.asarray(data_mat)
        # print(data_mat)
        final_data = np.zeros((num_rows//2,3))

        for i in range(num_rows//2):
            final_data[i,0] = data_mat[i,0]
            final_data[i,1] = data_mat[i,1]
            final_data[i,2]=  np.abs(data_mat[(i+num_rows//2),1] -final_data[i,1])
            # print(pos,mean,std)
            # final_data[i] = [pos,mean,sem]
        # print(final_data)
        return final_data

def get_valid_float(x):
    if not x == '':
        return float(x)
    return np.nan
def ReadCSVTime(filename):
    """
   function to read CSV data files:
       arguments : filename : path of the .csv file
   """
    data_mat = []
    with open(filename,encoding='utf-8-sig') as csvDataFile:
        # read file as csv file
        csvReader = csv.reader(csvDataFile)
        # loop over rows
        num_rows = 0
        for row in csvReader:
            # add cell [0] to list of dates
            if row != []:
                num_rows += 1
                # breakpoint()
                row_data = [float(x) for x in row]
                data_mat.append(row_data)
        data_mat = np.asarray(data_mat)
        return data_mat
def patterson_data_1g():
    data = ReadCSVFull(Patterson_2010_1g_file)
    return data
def patterson_data_4c():
    data = ReadCSVFull(Patterson_2010_4c_file)
    return data

def patterson_data_1b():
    data = ReadCSVFull(Patterson_2010_1b_file)
    return data
def tanaka_data_s3():
    psd_data = ReadCSVFull(Tanaka_2012_S3a_PSD_file)
    dend_data = ReadCSVFull(Tanaka_2012_S3a_Dend_file)
    return psd_data, dend_data

def graves_data_5d():
    dsred_data = ReadCSVTime(Graves_2021_dsred)
    glua1_data = ReadCSVTime(Graves_2021_glua1)
    return glua1_data,dsred_data

def get_tanaka_exo_data():
    exo_data = ReadCSVFull(Tanaka_2012_exo_data)
    return exo_data


# patterson_spine_sizes = patterson_data_1b()
# breakpoint()
# glua1_data,dsred_data = graves_data_5d()
# ng1d = glua1_data.T/glua1_data[:,0]
# mean_g1 = ng1d.mean(axis=1)
# tps = np.array([-10,-5,0,10,20,30,40,50,60])
# sem_g1 = ng1d.std(axis=1)/np.sqrt(glua1_data.shape[0])
# dsredd = dsred_data.T/dsred_data[:,0]
# mean_dsred = dsredd.mean(axis=1)
# sem_dsred = dsredd.std(axis=1)/np.sqrt(dsred_data.shape[0])
# # breakpoint()
# for i in range(glua1_data.shape[0]):
#     plt.plot(tps,ng1d[:,i],'k',alpha=0.2)
# plt.errorbar(tps,mean_g1,sem_g1,color="green",label="SEP-gluA1")
# # plt.errorbar(tps,mean_dsred,sem_dsred,color="red",label="dsRed")
# plt.hlines(y=1,xmin=tps[0],xmax=tps[-1])
# plt.legend()
# plt.xlabel("Time relative to induction (min)")
# plt.ylabel("Normalized signal (a.u.)")
# plt.xticks(tps)
# plt.show()


# exo_data = get_tanaka_exo_data()
# breakpoint()
def ReadCSVWillig(filename,head_rows = 0):
    """
   function to read CSV data files:
       arguments : filename : path of the .csv file
   """
    data_mat = []
    with open(filename) as csvDataFile:
        # read file as csv file
        csvReader = csv.reader(csvDataFile)
        # loop over rows
        num_rows = 0
        for rid,row in enumerate(csvReader):
            # add cell [0] to list of dates
            if row != [] and rid > head_rows-1:
                num_rows += 1
                row_data = [get_valid_float(x) for x in row[0:4]]
                data_mat.append(row_data)
        final_data = np.asarray(data_mat)
            # final_data[i] = [pos,mean,sem]
        # print(final_data)
        return final_data

def get_CF_2023_4_data():
    # data = {}
    # for fname in Clavet_Fournier_2023_4_fnames:
    #     file = fname+".csv"
    #     data[fname] = ReadCSVWillig(os.path.join(Clavet_Fournier_2023_4_folder,file),2)
    # breakpoint()
    control_timepoints_dict = {'Control_0min' : 0,
                       'Control+ 30 min' : 30,
                       'Control+ 1h' : 60,
                       'Control+ 2h' : 120, }
    cLTP_timepoints_dict = {'cLTP_0 min' : 0,
                       'cLTP+ 30 min' : 30,
                       'cLTP+ 60 min' : 60,
                       'cLTP + 2h' : 120}

    df = pd.read_excel(os.path.join(root_folder,'K_willig/TableS4_relto_Figure4.xlsx'), sheet_name=None)
    print(df.keys())
    control_df = pd.DataFrame()
    cLTP_df = pd.DataFrame()
    for ct_key in control_timepoints_dict.keys():
        temp_df = pd.DataFrame()
        temp_df['PSD_Area'] = df[ct_key].iloc[1:,0]
        temp_df['Morphology'] = df[ct_key].iloc[1:, 1]
        temp_df['GluA2_Area'] = df[ct_key].iloc[1:, 2]
        temp_df['#_Cluster'] = df[ct_key].iloc[1:, 3]
        temp_df['timepoint'] = control_timepoints_dict[ct_key]
        temp_df['condition'] = 'Control'
        control_df = pd.concat((control_df,temp_df)).dropna()
    for ct_key in cLTP_timepoints_dict.keys():
        temp_df = pd.DataFrame()
        temp_df['PSD_Area'] = df[ct_key].iloc[1:,0]
        temp_df['Morphology'] = df[ct_key].iloc[1:, 1]
        temp_df['GluA2_Area'] = df[ct_key].iloc[1:, 2]
        temp_df['#_Cluster'] = df[ct_key].iloc[1:, 3]
        temp_df['timepoint'] = cLTP_timepoints_dict[ct_key]
        temp_df['condition'] = 'cLTP'
        cLTP_df = pd.concat((cLTP_df,temp_df)).dropna()

    # dataframe2 = dataframe2.iloc[3:]
    return pd.concat((control_df,cLTP_df))
# CF_PSD = get_CF_2023_4_data()
# CF_control_psd = CF_PSD.loc[CF_PSD['condition'] == "Control"].dropna()
# CF_cLTP_psd = CF_PSD.loc[CF_PSD['condition'] == "cLTP"].dropna()
# CF_time_points = CF_PSD["timepoint"].unique()
# COI = "GluA2_Area"
#
# cLTP_means = CF_cLTP_psd.groupby('timepoint')[[COI]].mean().to_numpy()
# cLTP_means /= cLTP_means[0]
# cLTP_stds = CF_cLTP_psd.groupby('timepoint')[[COI]].std().to_numpy()
# cLTP_counts = CF_cLTP_psd.groupby('timepoint')[[COI]].count().to_numpy()
# cLTP_sems = cLTP_stds/cLTP_counts
#
# cont_means = CF_control_psd.groupby('timepoint')[[COI]].mean().to_numpy()
#     # breakpoint()
# cont_means /= cont_means[0]
# cont_stds = CF_control_psd.groupby('timepoint')[[COI]].std().to_numpy()
# cont_counts = CF_control_psd.groupby('timepoint')[[COI]].count().to_numpy()
# cont_sems = cont_stds / cont_counts
# fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
# print(ax.shape)
# cltp_sizes = CF_cLTP_psd.groupby('timepoint').GluA2_Area.apply(list).reset_index()[[COI]].to_numpy()
# control_sizes = CF_control_psd.groupby('timepoint').GluA2_Area.apply(list).reset_index()[[COI]].to_numpy()
# histtype='step'
# nbin=30
# for i in range(0,4):
#     plt.title(COI)
#     # ax[i//2,i%2].set_xscale('log')
#     ax[i//2,i%2].hist(np.log(cltp_sizes[i,0]),bins=nbin,label="ClTP",histtype=histtype,density=True,cumulative=True)
#     p1 = sp.stats.mstats.normaltest(np.log(cltp_sizes[i,0])).pvalue
#     ax[i // 2, i % 2].hist(np.log(control_sizes[i, 0]), bins=nbin, label="control",histtype=histtype,density=True,cumulative=True)
#     p2 = sp.stats.mstats.normaltest(np.log(control_sizes[i,0])).pvalue
#     print(p1,p2)
#     ax[i//2,i%2].set_title("time = {} mins".format(CF_time_points[i]))
#     ax[i//2,i%2].legend()
# plt.show()
# plt_widget = SNSPlottingWidget()
# save_it = 0
# y="GluA2_Area"
# hue="condition"
# x="timepoint"
# labs = ["Control","cLTP"]
# pairs = [
#     [(0,"Control"),(0,"cLTP")],
#     [(1,"Control"),(1,"cLTP")],
#     [(2,"Control"),(2,"cLTP")],
#     [(3,"Control"),(3,"cLTP")],
#     # [("Control",1),("cLTP",1)],
#     # [("Control",2),("cLTP",2)],
#     # [("Control",3),("cLTP",3)]
#
# ]
# #,("TTX","control"),("control","Untreated")]
# xlab = "Condition"
# ylab = r"$\frac{sGlua2}{iGluA2}$"
# title = ""
# stat = "GluA2_Area"
# fig_file = os.path.join(Clavet_Fournier_2023_4_folder,"{}_cLTP_temporal_profile".format(stat))
# stat_test = "Kruskal"
# fsize = 18
# order = ['Control','cLTP']
# states_order = ['Control','cLTP']
# subcat_order = [i for i in range(0,4)]
# subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
# states_palette = sns.color_palette("YlGnBu", n_colors=2)
# plt_widget.Swarmboxplotcombo(data=All_data,
#                              x=x,
#                              y=y,
#                              hue=hue,
#                              xlab=xlab,
#                              ylab=ylab,
#                              pairs=pairs,
#                              hue_order = states_order,
#                              order = subcat_order,
#                              labs=labs,
#                              title=title,
#                              color_pal=states_palette,
#                              stat_test=stat_test,
#                              xfsize=fsize,
#                              yfsize=1.2*fsize,
#                              fname=fig_file,
#                              save_it = save_it)
#
# fig, ax =  plt.subplots(figsize=(8,6),ncols=1,nrows=1)
# x_reg = All_data.dropna()
# sns.regplot(data=x_reg[x_reg["condition"]=="cLTP"],
#            x=x,
#            y=y)
# # sns.
# plt.show()