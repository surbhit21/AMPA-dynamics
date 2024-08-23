import matplotlib.pyplot as plt
# from matplotlib import rc
import matplotlib

from Utility import _2gaussian, _1gaussian

matplotlib.use("Qt5Agg")
import pandas as pd
from scipy.stats import ks_2samp, kruskal,spearmanr,pearsonr,ttest_rel,ttest_ind
import seaborn as sns
from SNSPlottingWidget import SNSPlottingWidget
from statannotations.Annotator import Annotator
from Utility import *
# from SpineEnrichment_fitting import FitModel
scale = 61.57/2048
# COLORS = {"spine":"#005f73","shaft":'#CA6702',"spine_s":"#005f73","spine_i":'#CA6702',"shaft_s":"#005f73","shaft_i":'#CA6702'}

my_pal1 = {"spine":COLORS_dict["soma"],"shaft":COLORS_dict["shaft"]}
my_pal2 = {"int": COLORS_dict["spine_i"],"int_getz": COLORS_dict["shaft_i"],"surf": COLORS_dict["spine_s"],"surf_getz": COLORS_dict["shaft_s"]}
my_pal3 = {"int": COLORS_dict["shaft_i"],"surf": COLORS_dict["shaft_s"]}
my_pal4 = {"Helm et al. 2021": COLORS_dict["spine_s"],"Getz et al. 2022": COLORS_dict["spine_i"]}




def CorrelationCalculationAndPlotting(data_to_show,localization,compartment,x_param,y_param,y_bg_param,f_scale=1):
    fig,ax = plt.subplots(figsize=(10,6),ncols=2,nrows=1)
    sns.set(font_scale = f_scale)
    sns.scatterplot(x=x_param,y=y_param,data=data_to_show,ax=ax[0],hue='dataset').set(title="Area vs {} intensity in {} rois".format(subunit,localization))
    sns.scatterplot(x=x_param,y=y_bg_param,data=data_to_show,ax=ax[1],hue='dataset').set(title="Area vs {} intensity in {} background rois".format(subunit,localization))
    sns.regplot(x=x_param,y=y_param,data=data_to_show[data_to_show["dataset"] ==compartment[0]],ax=ax[0])
    sns.regplot(x=x_param,y=y_param,data=data_to_show[data_to_show["dataset"] ==compartment[1]],ax=ax[0])
    sns.regplot(x=x_param,y=y_bg_param,data=data_to_show[data_to_show["dataset"] ==compartment[0]],ax=ax[1])
    sns.regplot(x=x_param,y=y_bg_param,data=data_to_show[data_to_show["dataset"] ==compartment[1]],ax=ax[1])
    corr1 = spearmanr(data_to_show[data_to_show['dataset']==compartment[0]][x_param],data_to_show[data_to_show['dataset']==compartment[0]][y_param])
    corr2 = spearmanr(data_to_show[data_to_show['dataset']==compartment[0]][x_param],data_to_show[data_to_show['dataset']==compartment[0]][y_bg_param])
    corr3 = spearmanr(data_to_show[data_to_show['dataset']==compartment[1]][x_param],data_to_show[data_to_show['dataset']==compartment[1]][y_param])
    corr4 = spearmanr(data_to_show[data_to_show['dataset']==compartment[1]][x_param],data_to_show[data_to_show['dataset']==compartment[1]][y_bg_param])
    print("Correlation between {} , {} in {} {} = ".format(x_param,y_param,compartment[0],localization),corr1)
    print("Correlation between {} , {} in {} {} = ".format(x_param,y_bg_param,compartment[0],localization),corr2)
    print("Correlation between {} , {} in {} {} = ".format(x_param,y_param,compartment[1],localization),corr3)
    print("Correlation between {} , {} in {} {} = ".format(x_param,y_bg_param,compartment[1],localization),corr4)

subunit = "GluA2"
p_type="SUM"
stat_test = "Kruskal"
condition = "control" # or "
# TTX" or "Untreated"
    # {}

G2std = GluA2StoD(spine_cell_folder[subunit][condition],p_type,spine_channel_names[subunit])
int_df,surf_df,gfp_df,\
int_dend_data,surf_dend_data,gfp_dend_data,\
int_spine_line_data,surf_spine_line_data,gfp_spine_line_data = G2std.LoadData(
    spine_num_cells[subunit][condition],
    spine_exclude_cells[subunit][condition])
concat_data = pd.concat([int_df.assign(dataset='internal'),surf_df.assign(dataset='surface')])
# print(concat_data.shape)
# breakpoint()

# concat_data
# int_dend_data
# surf_df#[surf_df.index == "cell_9_sp_10"]
# breakpoint()
plt_widget = SNSPlottingWidget()
op_folder = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/{}/Figures/Spine_to_dend/".format(subunit)
plt_widget.CreateFolderRecursive(op_folder)
print(pd. __version__)
ax_label = 1 #plot axis labels or not (=1 for yes)
fsize = 28
save_it = 1
w_or_wo_ax_label = ["", "_with_ax_label"] #for chaning the filename suffix

## setting for all hist plots
cumm = False
stat_to_plot = "probability"
log_scale= True
legend = True
element="poly"
legend = True
alpha = 1
common_norm = False
fill = False
# sns.set(font_scale = 2.5)
"""
## Sanity check to see if the data is loaded correctly
# 1. Area of ROI and background ROI should be the same
roi_area_ratio = pd.DataFrame()
roi_area_ratio['spine'] = int_df["sp_area"]/int_df["sp_bg_area"]
roi_area_ratio['shaft'] = int_df["dend_area"]/int_df["dend_bg_area"]
fig,ax = plt.subplots(figsize=(12,8),ncols=1,nrows=1)
# ratio_of_ratio = ratio_SE.eval("int / surf").rename("Ratio of Int/surf ratio")
x1 = np.ones(roi_area_ratio['spine'].size)
sns.scatterplot(data=roi_area_ratio,ax=ax).set(xticklabels=[])
y1 = x1
sns.lineplot(x1,y1,color='r',label="y=1")
# roi_area_ratio[roi_area_ratio['shaft']<1]
plt_widget.SaveFigures(os.path.join(op_folder,"are_ratio"))
plt.show()
"""

"""
spine_areas = int_df["sp_area"]
dend_areas = int_df["dend_area"]
fig,ax1 = plt.subplots()

hue = "condition"
loc = "surface"

sp1 = sns.histplot(int_df["sp_area"],
                   bins=10,
                   cumulative=True,
                   alpha=alpha,
                   element=element,
                   stat=stat_to_plot,
                   common_norm=False,
                   fill=fill,
                   legend=legend,
                   color=COLORS_dict["spine"],
                   ax=ax1)
sp2 = sns.histplot(int_df["dend_area"],
                   bins=10,
                   cumulative=True,
                   alpha=alpha,
                   element=element,
                   stat=stat_to_plot,
                   common_norm=False,
                   fill=fill,
                   legend=legend,
                   color=COLORS_dict["shaft"],
                   ax=ax1)
if ax_label==1:
        # breakpoint()
    ax1.set_xlabel(r"Spine area ($\mu m^2$)",fontsize=fsize)
    ax1.set_ylabel(stat_to_plot,fontsize=fsize)
    ax1.set_title("spine and shaft area histogram")
ax1.legend(["spine","shaft"])
plt.show()
"""

# breakpoint()
"""
## ratio of spine to dendrite area
s_2_d_area_ratio = pd.DataFrame()
s_2_d_area_ratio['area_ratio'] = int_df["sp_area"]/int_df["dend_area"]
print(s_2_d_area_ratio)
s_2_d_area_ratio.iteritems = s_2_d_area_ratio.items
sns.boxplot(data=s_2_d_area_ratio)
# outliers = [y for stat in boxplot_stats(s_2_d_area_ratio) for y in stat['fliers']]
# outlier_area_df = s_2_d_area_ratio[s_2_d_area_ratio['area_ratio'].isin(outliers)]
# print(outlier_area_df)
plt.show()

## getting number of pixels in each roi

num_pixels_data = pd.DataFrame()
num_pixels_data['spine'] = surf_df['sp_area']
num_pixels_data['shaft'] = surf_df['dend_area']
# dend_pixels = surf_df['dend_area']/0.0009
fig,ax = plt.subplots(figsize=(12,8),ncols=1,nrows=1)
sns.scatterplot(x="spine",y="shaft",data=num_pixels_data,ax=ax).set(title="Area of spine vs shaft")
sns.regplot(x="spine",y="shaft",data=num_pixels_data,ax=ax)
x1 = np.arange(0,num_pixels_data['spine'].max(),0.001)
y1 = x1
# sns.lineplot(x1,y1,color="r",label="x=y")
# sp_pixels,dend_pixels
plt.show()

## Correlation of area and GluA2 intensity in spines vs background

stat1 = "mean"
stat2 = "area"
compartments = ["surface","internal"]
localization = ["sp","dend"]
y_param = "{}_{}".format(localization[0],stat1)
x_param = "{}_{}".format(localization[0],stat2)
y_bg_param = "{}_bg_{}".format(localization[0],stat1)
dts = concat_data
CorrelationCalculationAndPlotting(dts,"spine",compartments,x_param,y_param,y_bg_param,f_scale=1)
y_param = "{}_{}".format(localization[1],stat1)
x_param = "{}_{}".format(localization[1],stat2)
y_bg_param = "{}_bg_{}".format(localization[1],stat1)
CorrelationCalculationAndPlotting(dts,"shaft",compartments,x_param,y_param,y_bg_param,f_scale=1)


## Histogram of mean intensity in spines and spine background ROIs to see if spines can be distinguished from background
num_bins = 10
print(condition)
fig,ax = plt.subplots(figsize=(12,8),ncols=2,nrows=1)
ax[0].set_title("spine internal")
ax[1].set_title("spine surface")
sns.histplot(int_df["sp_bg_mean"],bins=num_bins, alpha=0.8,
    stat="count", legend=True,color=COLORS_dict["spine_i"],ax=ax[0],fill=True)
sns.histplot(int_df["sp_mean"],bins=num_bins,  alpha=0.2,
    stat="count",legend=True,color=COLORS_dict["spine_i"],ax=ax[0],fill=False).set( xlabel=r"Integral density ",ylabel=r"Count")
sns.histplot(surf_df["sp_bg_mean"],bins=num_bins, alpha=0.8,
    stat="count",legend=True,color=COLORS_dict["spine_s"],ax=ax[1],fill=True)
sns.histplot(surf_df["sp_mean"],bins=num_bins, alpha=0.2,
    stat="count",legend=True,color=COLORS_dict["spine_s"],ax=ax[1],fill=False).set( xlabel=r"Integral density ",ylabel=r"Count")
ax[0].legend(["background","signal"],bbox_to_anchor = (1, 1))
ax[1].legend(["background","signal"],bbox_to_anchor = (1, 1))
# ax.errorbar(np.arange(0,17,1),int_df['sp_mean'],int_df['sp_stddev'])
plt.show()

## Histogram of mean intensity in shaft and shaft background to see if shafts can be distinguished from background

fig,ax = plt.subplots(figsize=(12,8),ncols=2,nrows=1)
ax[0].set_title("shaft internal")
ax[1].set_title("shaft surface")
sns.histplot(int_df["dend_bg_mean"],bins=num_bins, alpha=0.8,
    stat="count", legend=True,color=COLORS_dict["shaft_i"],ax=ax[0])
sns.histplot(int_df["dend_mean"],bins=num_bins,  alpha=0.2,
    stat="count",legend=True,color=COLORS_dict["shaft_i"],ax=ax[0],fill=False).set( xlabel=r"Integral density ",ylabel=r"Count")
sns.histplot(surf_df["dend_bg_mean"],bins=num_bins, alpha=0.8,
    stat="count",legend=True,color=COLORS_dict["shaft_s"],ax=ax[1])
sns.histplot(surf_df["dend_mean"],bins=num_bins, alpha=0.2,
    stat="count",legend=True,color=COLORS_dict["shaft_s"],ax=ax[1],fill=False).set( xlabel=r"Integral density ",ylabel=r"Count")
ax[0].legend(["background","signal"],bbox_to_anchor = (1, 1))
ax[1].legend(["background","signal"],bbox_to_anchor = (1, 1))
plt.show()

"""


def GetSEvals(stat,s_df,i_df,norm_by_area):
    ratio_SE = pd.DataFrame()
    if norm_by_area == True:
        ratio_SE["surf"] = ((s_df['sp_{}'.format(stat)] - s_df['sp_bg_{}'.format(stat)]) / s_df["sp_area"]) / (
                (s_df['dend_{}'.format(stat)] - s_df['dend_bg_{}'.format(stat)]) / s_df["dend_area"])
        ratio_SE["int"] = ((i_df['sp_{}'.format(stat)] - i_df['sp_bg_{}'.format(stat)]) / s_df["sp_area"]) / (
                (i_df['dend_{}'.format(stat)] - i_df['dend_bg_{}'.format(stat)]) / s_df["dend_area"])
    else:
        ratio_SE["surf"] = ((s_df['sp_{}'.format(stat)] - s_df['sp_bg_{}'.format(stat)])) / (
            (s_df['dend_{}'.format(stat)] - s_df['dend_bg_{}'.format(stat)]))
        ratio_SE["int"] = ((i_df['sp_{}'.format(stat)] - i_df['sp_bg_{}'.format(stat)])) / (
            (i_df['dend_{}'.format(stat)] - i_df['dend_bg_{}'.format(stat)]))
    ratio_SE["id"] = s_df.index
    return ratio_SE

def calculateSE(stat,s_df,i_df,x,y,pairs,labs,order,hue_order,fname,x_pos = [0,1],norm_by_area = True):
    ratio_SE = GetSEvals(stat,s_df,i_df,norm_by_area)
    y_pos = ratio_SE[["surf", "int"]].to_numpy()
    ratio_SE_melted = ratio_SE.melt(id_vars=["id"], var_name="localization", value_name="SE")
    plt_widget.SwarmboxLineplotcombo(data=ratio_SE_melted,
                                     x=x,
                                     y=y,
                                     hue=[],
                                     pairs=pairs,
                                     xlab="Compartment",
                                     ylab="Synaptic enrichment",
                                     labs=labs,
                                     title="",
                                     order=order,
                                     hue_order=hue_order,
                                     color_pal=my_pal3,
                                     stat_test=stat_test,
                                     xfsize=fsize,
                                     yfsize=fsize,
                                     fname=fig_file,
                                     save_it=save_it,
                                     ax_lab=ax_label,
                                     x_pos=x_pos,
                                     line_y_data=y_pos,
                                     )
## spine to dendrite ratio using the method from Helm et al. 2021
# """

"""
plotting parameters are set
"""
x="localization"
y = "SE"
pairs = [("surf","int")]
labs = ["Spine","Shaft"]
order = ["surf","int"]
hue = order
hue_order = order
# breakpoint()
# breakpoint()
x_pos = [0,1]


"""
SE based on intdent, rawinteden and mean  
"""
stat = "intden"
fig_file = os.path.join(op_folder,"spine_to_dend_{0}_ratio_{1}".format(stat,w_or_wo_ax_label[ax_label]))
calculateSE(stat,surf_df,int_df,x,y,pairs,labs,order,hue_order,fig_file,norm_by_area=False)

stat = "rawintden"
fig_file = os.path.join(op_folder,"spine_to_dend_{0}_ratio_{1}".format(stat,w_or_wo_ax_label[ax_label]))
calculateSE(stat,surf_df,int_df,x,y,pairs,labs,order,hue_order,fig_file,norm_by_area=False)

stat = "mean"
fig_file = os.path.join(op_folder,"spine_to_dend_{0}_ratio_{1}".format(stat,w_or_wo_ax_label[ax_label]))
calculateSE(stat,surf_df,int_df,x,y,pairs,labs,order,hue_order,fig_file,norm_by_area=False)


"""
plotting cummulative histogram of SE values
"""

stat = "intden"
ratio_SE = GetSEvals(stat,surf_df,int_df,norm_by_area=True)
ratio_SE_melted = ratio_SE.melt(id_vars=["id"],var_name="localization",value_name="SE")
fig,ax1 = plt.subplots()
hp1 = sns.histplot(ratio_SE_melted,
             x=y,
             hue=x,
             hue_order = ["surf","int"],
             bins=10,
             cumulative=cumm,
             alpha=alpha,
             element=element,
             stat=stat_to_plot,
             common_norm=False,
             fill=fill,
             legend=legend,
             color=COLORS_dict["spine"],
             ax=ax1).set(label="S-{}".format(subunit))#set()
# hp2 = sns.histplot(ratio_SE["int"],
#              bins=10,
#              cumulative=True,
#              alpha=alpha,
#              element=element,
#              stat=stat_to_plot,
#              common_norm=False,
#              fill=fill,
#              legend=legend,
#              color=COLORS_dict["spine_i"],
#              ax=ax1).set(label="I-{}".format(subunit))#set(ylabel="Cummulative frequency",xlabel="Synaptic enrichment")
if ax_label==1:
    ax1.set_xlabel("Synaptic enrichment", fontsize=fsize)
    ax1.set_ylabel(ylabel="Cummulative frequency", fontsize=fsize)
    ax1.set_title( "Synaptic enrichment histogram")
    legend = ax1.get_legend()
    handles = legend.legendHandles
    legend.remove()
    ax1.legend(handles, ["S-{}".format(subunit), "i-{}".format(subunit)], title='',frameon=False,fontsize=fsize)
# ax1.legend()
ax1.set_ylim([0,1.1])
ax1.tick_params(axis='both', which='major', labelsize=fsize)
ax1.grid(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# ax1.legend()
plt_widget.SaveFigures(os.path.join(op_folder,"synaptic_enrichment_histogram_{0}".format(w_or_wo_ax_label[ax_label])))
plt.show()


# adding another column named binned with bin number
bin_size = 10
bins = np.arange(0,61,bin_size)
labels = ((bins/10)+1)[:-1]
labels = labels.astype(int)
# bins,labels
binned_spines = pd.DataFrame()
concat_data["binned"] = pd.cut(concat_data["DFO"],bins=bins,labels=labels)
surf_df["binned"] = pd.cut(surf_df["DFO"],bins=bins,labels=labels)
int_df["binned"] = pd.cut(int_df["DFO"],bins=bins,labels=labels)
gfp_df["binned"] = pd.cut(gfp_df["DFO"],bins=bins,labels=labels)

fig,ax5 = plt.subplots(figsize=(8,6),nrows=1,ncols=1)
# ax5.grid(False)
ax5.tick_params(axis='both', which='major', labelsize=fsize)
ax5.spines['right'].set_visible(False)
ax5.spines['top'].set_visible(False)
# ax5[1].spines['right'].set_visible(False)
# ax5[1].spines['top'].set_visible(False)
# ax5.labels

stat = "intden"
s2s_r_dist = GetSEvals(stat,surf_df,int_df,norm_by_area=False)
# s2s_r_dist["int"]  = ((int_df['sp_{}'.format(stat)]-int_df['sp_bg_{}'.format(stat)]))/((int_df['dend_{}'.format(stat)] - int_df['dend_bg_{}'.format(stat)]))
# s2s_r_dist["surf"] = ((surf_df['sp_{}'.format(stat)]-surf_df['sp_bg_{}'.format(stat)]))/((surf_df['dend_{}'.format(stat)] - surf_df['dend_bg_{}'.format(stat)]))
s2s_r_dist["binned"] = surf_df["binned"]
s2s_r_dist["DFO"] = surf_df["DFO"]
y_tics = [str(int(10*i)) for i in labels]
# s2s_r_dist.melt()


# sns.stripplot(data=s2s_r_dist,y="surf",x="binned",ax=ax5[0]).set(xlabel=r"Dendritic distance ($\mu m$)",ylabel=r"ratio $\frac{spine}{shaft}$")
# sns.swarmplot(data=s2s_r_dist,y="int",x="binned",ax=ax[1]).set(xlabel=r"Dendritic distance ($\mu m$)",ylabel=r"ratio $\frac{spine}{shaft}$")
# sns.swarmplot(data=mean_s2s_r_dist,y="surf",x="binned",ax=ax[0]).set(xlabel=r"Dendritic distance ($\mu m$)",ylabel=r"ratio $\frac{spine}{shaft}$")
# sns.swarmplot(data=mean_s2s_r_dist,y="int",x="binned",ax=ax[1]).set(xlabel=r"Dendritic distance ($\mu m$)",ylabel=r"ratio $\frac{spine}{shaft}$")
sns.regplot(data=s2s_r_dist,y="surf",x="DFO",ax=ax5,scatter_kws={"color":my_pal3["surf"],"alpha" : 0.5},line_kws={"color": my_pal3["surf"],"alpha":1.0}).set(xlabel="",ylabel="")
# sns.regplot(data=s2s_r_dist,y="int",x="DFO",ax=ax5[1], scatter_kws={"color":my_pal3["int"],"alpha" : 1.0},line_kws={"color": my_pal3["int"]}).set(ylabel="",xlabel="")#.set(xlabel=r"Dendritic distance ($\mu m$)",ylabel=r"ratio $\frac{spine}{shaft}$")
if ax_label==1:
   ax5.set_xlabel(r"Dendritic distance ($\mu m$)",fontsize=fsize)
   ax5.set_ylabel(r"Synaptic enrichment",fontsize=fsize)
pea_ror_s = spearmanr(s2s_r_dist["surf"],s2s_r_dist["DFO"])
pea_ror_i = spearmanr(s2s_r_dist["int"],s2s_r_dist["DFO"])
# print(pea_ror)
ax5.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(pea_ror_s[0],pea_ror_s[1]),color=my_pal3["surf"],fontsize=28,transform=ax5.transAxes)
ax5.set_xlim([0,60])
# ax5[1].text(.05, .8, 'r={:.2f}, p={:.2g}'.format(pea_ror_i[0],pea_ror_i[1]),transform=ax5[1].transAxes)

# ax[0].set_xlim([0,7.5])
# print(y_tics)
titles = ["Surface","Internal"]
# for adx,axes in enumerate(ax5):
#     # axes.set_xticklabels(y_tics)
#     if ax_label==1:
#         axes.set_xlabel(r"Dendritic distance [$\mu m$]")
#         axes.set_ylabel(r"ratio $\frac{spine}{shaft}$")
#         axes.set_title(titles[adx])
#     axes.set_xlim([0,60])
plt_widget.SaveFigures(os.path.join(op_folder, "surf_spine_to_dend_{0}_ratio_distribution_{1}".format(stat,w_or_wo_ax_label[ax_label])))
plt.show()
# s2s_r_dist

## fitting mode to spine to shaft ratio dendritic distribution

# mean_s2s_r_dist = s2s_r_dist.groupby(['binned'],as_index=False).mean()
# std_s2s_r_dist = s2s_r_dist.groupby(['binned'],as_index=False).std()
# x_arr = np.asarray(mean_s2s_r_dist["binned"])*10
# se_arr = np.asarray(mean_s2s_r_dist["surf"])
# std_se_arr = np.asarray(std_s2s_r_dist["surf"])
# # breakpoint()
# x1,SE_dist,SE_chi2,paras,mini,out2 = FitModel(x_arr,se_arr,std_se_arr,0.4,np.Inf,60)
# # fig,ax = plt.subplots()
# # sns.swarmplot(data=s2s_r_dist,y="surf",x="binned",ax=ax).set(xlabel=r"Distance from Soma ($\mu m$)",ylabel=r"ratop $\frac{spine}{shaft}$",label="data")
# sns.lineplot(x1,SE_dist,ax=ax5[0]).set(label=r"$\frac{P_{spine}}{P_s}$model-fit, $\chi^2$=%0.2f" %(SE_chi2))
# plt.show()
# breakpoint()

"""
distribution of inteden as a function of distance for spine and dend seperately
"""
fig,ax6 = plt.subplots(figsize=(8,6),nrows=1,ncols=1)
ax6 = [ax6]
for ax in ax6:
    ax.tick_params(axis='both', which='major', labelsize=fsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

surf_df_corrected = pd.DataFrame()
int_df_corrected = pd.DataFrame()
stat = "intden"
surf_df_corrected["spine"] = ((surf_df['sp_{}'.format(stat)]-surf_df['sp_bg_{}'.format(stat)]))
surf_df_corrected["shaft"] = ((surf_df['dend_{}'.format(stat)] - surf_df['dend_bg_{}'.format(stat)]))
surf_df_corrected["binned"] = surf_df["binned"]
surf_df_corrected["DFO"] = surf_df["DFO"]
int_df_corrected["spine"] = ((int_df['sp_{}'.format(stat)]-int_df['sp_bg_{}'.format(stat)]))
int_df_corrected["shaft"] = ((int_df['dend_{}'.format(stat)] - int_df['dend_bg_{}'.format(stat)]))
int_df_corrected["binned"] = int_df["binned"]
int_df_corrected["DFO"] = int_df["DFO"]
y_tics = [str(int(10*i)) for i in labels]
sns.regplot(data=surf_df_corrected,y="spine",x="DFO",ax=ax6[0],marker = "o",label="S-spine",scatter_kws={"color": my_pal3["surf"],"alpha" : 0.5},line_kws={"color": my_pal3["surf"],'alpha':1.0})
sns.regplot(data=surf_df_corrected,y="shaft",x="DFO",ax=ax6[0],marker = "v",label="S-shaft", scatter_kws={"color":my_pal3["int"],"alpha" : 0.5},line_kws={"color": my_pal3["int"],'alpha':1})
# sns.regplot(data=int_df_corrected,y="spine",x="DFO",ax=ax6[1],label="I-spine",scatter_kws={"color":my_pal2["surf"],"alpha" : 1.0},line_kws={"color": my_pal2["surf"]})
# sns.regplot(data=int_df_corrected,y="shaft",x="DFO",ax=ax6[1],label="I-shaft", scatter_kws={"color":my_pal2["int"],"alpha" : 1.0},line_kws={"color": my_pal2["int"]})
# ax6[0].legend(title='localization', loc='upper right', labels=['spine',None,None, 'shaft'])
# ax6[1].legend(title='localization', loc='upper right', labels=['spine','','', 'shaft'])

pea_ror_s = spearmanr(surf_df_corrected["spine"],s2s_r_dist["DFO"])
pea_ror_i = spearmanr(surf_df_corrected["shaft"],s2s_r_dist["DFO"])

ax6[0].text(.05, .8, 'r={:.2f},\np={:.2g}'.format(pea_ror_s[0],pea_ror_s[1]),color=my_pal3["surf"],fontsize=28,transform=ax6[0].transAxes)
ax6[0].text(.05, .8, 'r={:.2f},\np={:.2g}'.format(pea_ror_i[0],pea_ror_i[1]),color=my_pal3["int"],fontsize=28,transform=ax6[0].transAxes)
# ax6[0].set_xlim([0,60])
print(pea_ror_s,pea_ror_i)
titles = ["Surface","Internal"]
for adx,axes in enumerate(ax6):
    # axes.set_xticklabels(y_tics)
    if ax_label==1:
        axes.set_xlabel(r"Dendritic distance ($\mu m$)",fontsize=fsize)
        axes.set_ylabel(r"Integrated density",fontsize=fsize)
        axes.set_title(titles[adx])
    axes.set_xlim([0,60])
    axes.legend(frameon=False)
plt_widget.SaveFigures(os.path.join(op_folder,"spine_and_dend_{0}_intensity_distribution_{1}".format(stat,w_or_wo_ax_label[ax_label])))
plt.show()




fig,ax6 = plt.subplots(figsize=(8,6),nrows=1,ncols=1)
ax6.tick_params(axis='both', which='major', labelsize=fsize)
ax6.spines['right'].set_visible(False)
ax6.spines['top'].set_visible(False)

sp_area_dist = pd.DataFrame()
sp_area_dist["area"] = surf_df['sp_area']
sp_area_dist["dend area"] = surf_df['dend_area']
sp_area_dist["binned"] = surf_df["binned"]
sp_area_dist["DFO"] = surf_df["DFO"]
y_tics = [str(int(10*i)) for i in labels]
sns.regplot(data=sp_area_dist,y="area",x="DFO",ax=ax6,scatter_kws={"color":my_pal3["surf"],"alpha" : 0.5},line_kws={"color": my_pal3["surf"],"alpha":1.})
# sns.regplot(data=sp_area_dist,y="dend area",x="DFO",ax=ax6,scatter_kws={"color":my_pal3["int"],"alpha" : 0.5},line_kws={"color": my_pal3["int"],"alpha":1.})

# sns.regplot(data=sp_area_dist,y="dend area",x="DFO",ax=ax6[1],scatter_kws={"color":my_pal3["int"],"alpha" : 1.0},line_kws={"color": my_pal3["int"]})#.set(xlabel=r"Dendritic distance ($\mu m$)",ylabel=r"ratio $\frac{spine}{shaft}$")
pea_area_s = spearmanr(sp_area_dist["area"],sp_area_dist["DFO"])
# pea_area_d = spearmanr(sp_area_dist["dend area"],sp_area_dist["DFO"])
# pea_area_d = spearmanr(sp_area_dist["dend area"],sp_area_dist["DFO"])

print(pea_area_s)
ax6.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(pea_area_s[0],pea_area_s[1]),transform=ax6.transAxes)
# ax6.text(.05, .9, 'r={:.2f}, p={:.2g}'.format(pea_area_d[0],pea_area_d[1]),transform=ax6.transAxes)
# print(y_tics)

## Spine areas as a function of distance
comparment = ["Spine"]
titles = ["Distribution of spine size","Distribution of shaft size"]
for adx,axes in enumerate([ax6]):
    if ax_label == 1:
        axes.set_xlabel(r"Dendritic distance ($\mu m$)",fontsize=fsize)
        axes.set_ylabel(r"{} area ($\mu m^2)$".format(comparment[adx]),fontsize=fsize)
        axes.set_title(titles[adx])
        axes.set_xlim([0,60])

plt_widget.SaveFigures(os.path.join(op_folder, "spine_size_distribution_{0}".format(w_or_wo_ax_label[ax_label])))
plt.show()
# s2s_r_dist

"""
## Ratio of ratios

fig,ax = plt.subplots(figsize=(12,8),ncols=1,nrows=1)
ax.tick_params(axis='both', which='major', labelsize=fsize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# ratio_of_ratio = ratio_SE.eval("int / surf").rename("Ratio of Int/surf ratio")
sns.regplot(x="int",y="surf",data=ratio_SE,ax=ax).set(title = r"$\frac{spine}{shaft}$ in surface vs intrnal {}".format(subunit))
x1 = np.arange(0,ratio_SE.max().max(),0.1)
y1 = x1
sns.lineplot(x1,y1,color='r',label="x=y")
pea_ror = spearmanr(ratio_SE["int"],ratio_SE["surf"])
print(pea_ror)
ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(pea_ror[0],pea_ror[1]),transform=ax.transAxes)
# ratio_of_ratio.mean()
plt_widget.SaveFigures(os.path.join(op_folder,"ratio_of_ratio_{0}".format(w_or_wo_ax_label[ax_label])))
plt.show()
"""
ratio_s_i = pd.DataFrame()
stat = "mean"
# breakpoint()
ratio_s_i["Spine"]  = ((surf_df['sp_{}'.format(stat)]-surf_df['sp_bg_{}'.format(stat)]))/((int_df['sp_{}'.format(stat)] - int_df['sp_bg_{}'.format(stat)]))
ratio_s_i["Shaft"] = ((surf_df['dend_{}'.format(stat)]-surf_df['dend_bg_{}'.format(stat)]))/((int_df['dend_{}'.format(stat)] - int_df['dend_bg_{}'.format(stat)]))
ratio_s_i["id"] = surf_df.index
fig_file = os.path.join(op_folder,"surf_to_int_{0}_in_spine_vs_shaft_{1}".format(stat,w_or_wo_ax_label[ax_label]))
# breakpoint()
y_pos = ratio_s_i[["spine","shaft"]].to_numpy()

ratio_s_i_melted = ratio_s_i.melt(id_vars=["id"],
                             var_name="localization",
                             value_name="ratio")
breakpoint()
pairs = [("spine","shaft")]
x="localization"
y = "ratio"
labs = []
xlab = "Compartment"
s_glua = "s{}".format(subunit)
i_glua = "i{}".format(subunit)
ylab= r"$\frac{"+s_glua+"}{"+i_glua+"}$"
plt_widget.SwarmboxLineplotcombo(data=ratio_s_i_melted,
                             x=x,
                             y=y,
                             hue=None,
                             xlab=xlab,
                             ylab=ylab,
                             title="",
                             labs=labs,
                             pairs=pairs,
                             color_pal=my_pal1,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=save_it,
                             ax_lab=ax_label,
                             x_pos=x_pos,
                             line_y_data=y_pos
                            )
# """
# fig,ax = plt.subplots(nrows=1,ncols=1)
# ax.tick_params(axis='both', which='major', labelsize=fsize)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # sns.set(font_scale = 2)
# sns.boxplot(data=ratio_s_i,notch=True,palette = my_pal1,ax=ax)
# if ax_label==1:
#     ax.set_xlabel("compartment",fontsize=fsize)
#     ax.set_ylabel(r"ratio $\frac{surface}{internal}$",fontsize=fsize)
#     ax.set_title("suface to internal {} intensity ratio".format(stat),fontsize=fsize)
# annotator = Annotator(ax, data=ratio_s_i.melt(),pairs=pairs,x="variable",y="value")
# annotator.configure(test=stat_test, fontsize = 38)
# annotator.apply_and_annotate()
# # print(ratio_s_i.melt())
# # outliers = [y for stat in boxplot_stats(ratio_s_i['shaft']) for y in stat['fliers']]
# # print(outliers)
# # outlier_s_i_df = ratio_s_i[ratio_s_i['shaft'].isin(outliers)]
# print(ttest_ind(ratio_s_i["shaft"] ,ratio_s_i["spine"] ))
# # outlier_s_i_df
# # concat_data[concat_data.index == "cell_5_sp_5"]
# print(ratio_s_i["spine"].mean(),ratio_s_i["shaft"].mean(),ratio_s_i["spine"].std(),ratio_s_i["shaft"].std())
# plt_widget.SaveFigures(os.path.join(op_folder,"surf_to_int_{}_in_spine_vs_shaft".format(stat)))
# plt.show()



## Synaptic enrichment as per defined by Getz et al. 2022 Sci. Adv.

plotting_spines = ["cell_1_sp_4","cell_6_sp_5"]
def SE_Getz(cell_id,data,final_data,comp,c="b"):
    num_spine = len(data[cell_id].keys())

    df_min = 0
    df_max = 5000
    sigma_max = 1
    all_fits_params = []
    np.random.seed(2023)
    op_folder2 = os.path.join(op_folder, "getz_fits")
    plt_widget.CreateFolderRecursive((op_folder2))
    # print(cell_id,num_spine)
    for i in range(1,num_spine+1):
        x_S = data[cell_id]["sp_{}".format(i)][:,0]
        y_S = data[cell_id]["sp_{}".format(i)][:,1]
        # print(y_S)
        y_S /= y_S.max()
        mid = x_S[-1]/2
        # print(x_S[-1])
        min_bounds = (0, y_S.min(),   0,  y_S.min(),  mid,   0)
        max_bounds = (1, y_S.max(), 0.5, y_S.max(), x_S[-1], 0.5)

        p_init = []
        for ps in range(len(min_bounds)):
            p_init.append(np.random.uniform(min_bounds[ps],max_bounds[ps]))
        popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, x_S, y_S, p0=p_init,bounds=(min_bounds,max_bounds),maxfev = 30000)
        sp = "sp_{}".format(i)
        spid = "_".join([cell_id,sp])
        if spid in plotting_spines:
        # print(spid)
            pars_1 = popt_2gauss[0:3]
            pars_2 = popt_2gauss[3:6]
            gauss_peak_1 = _1gaussian(x_S, *pars_1)
            gauss_peak_2 = _1gaussian(x_S, *pars_2)
            fig, ax = plt.subplots(figsize=(8, 6 ), ncols=1, nrows=1)
            ax.grid(False)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.plot(x_S, y_S, color=c, label="data")
            ax.plot(x_S, gauss_peak_1, ".r-", label="spine-fit",linewidth=0.5,alpha=0.5,markersize=3)
            ax.plot(x_S, gauss_peak_2, ".b-", label="shaft-fit",linewidth=0.5,alpha=0.5,markersize=3)
            ax.set_xlabel(r"Distance ($\mu$m)",fontsize=fsize-1)
            ax.set_ylabel("{}{} Normalized Intensity".format(comp,subunit),fontsize=fsize-1)
            ax.legend()
            fname = os.path.join(op_folder2, "getz_fit_{}_{}".format(comp, spid))
            plt_widget.SaveFigures(fname)
            plt.close()
        final_data["{}_spine".format(comp)].append(spid)
        final_data["{}_sp_amp".format(comp)].append(popt_2gauss[0])
        final_data["{}_sp_cen".format(comp)].append(popt_2gauss[1])
        final_data["{}_sp_sigma".format(comp)].append(popt_2gauss[2])
        final_data["{}_dend_amp".format(comp)].append(popt_2gauss[3])
        final_data["{}_dend_cen".format(comp)].append(popt_2gauss[4])
        final_data["{}_dend_sigma".format(comp)].append(popt_2gauss[5])
        # final_data["compartment".format(comp)].append(comp)
        # print(popt_2gauss)
        # all_fits_params.append(popt_2gauss)


    return final_data

#
getz_fitss = {}
getz_fitss["s_spine"] = []
getz_fitss["s_sp_amp"] = []
getz_fitss["s_sp_cen"] = []
getz_fitss["s_sp_sigma"] = []
getz_fitss["s_dend_amp"] = []
getz_fitss["s_dend_cen"] = []
getz_fitss["s_dend_sigma"] = []
getz_fitss["i_spine"] = []
getz_fitss["i_sp_amp"] = []
getz_fitss["i_sp_cen"] = []
getz_fitss["i_sp_sigma"] = []
getz_fitss["i_dend_amp"] = []
getz_fitss["i_dend_cen"] = []
getz_fitss["i_dend_sigma"] = []
for i in range(1,spine_num_cells[subunit][condition]+1):
    if i not in spine_exclude_cells[subunit][condition]:
        cell_toa = "cell_{}".format(i)
        getz_fitss = SE_Getz(cell_toa,surf_spine_line_data,getz_fitss,"s",c="#f003f0ff")
        # fitted_params_surf.append(f_p_s)
for i in range(1,spine_num_cells[subunit][condition]+1):
    if i not in spine_exclude_cells[subunit][condition]:
        cell_toa = "cell_{}".format(i)
        getz_fitss = SE_Getz(cell_toa,int_spine_line_data,getz_fitss,"i",c="#00d2d2ff")

getz_data = pd.DataFrame.from_dict(getz_fitss)
getz_data = getz_data.rename(columns={'s_spine':'spine'})
getz_data = getz_data.drop(columns=["i_spine"])
getz_data = getz_data.set_index('spine')

ratio_getz = pd.DataFrame()
ratio_getz["surf"] = getz_data["s_sp_amp"]/getz_data["s_dend_amp"]
ratio_getz["int"] = getz_data["i_sp_amp"]/getz_data["i_dend_amp"]
ratio_melted_getz = ratio_getz.melt(var_name="localization",value_name="SE")
# ratio_getz["id"] = getz_data.index
ratio_melted_getz["Ref"] = "Getz et al. 2022"

stat = "mean"
ratio_helm = pd.DataFrame()
ratio_helm["surf"] = ((surf_df['sp_{}'.format(stat)]-surf_df['sp_bg_{}'.format(stat)]))/((surf_df['dend_{}'.format(stat)] - surf_df['dend_bg_{}'.format(stat)]))
ratio_helm["int"]  = ((int_df['sp_{}'.format(stat)]-int_df['sp_bg_{}'.format(stat)]))/((int_df['dend_{}'.format(stat)] - int_df['dend_bg_{}'.format(stat)]))
ratio_melted_helm = ratio_helm.melt(var_name="localization",value_name="SE")
# ratio_helm["id"] = surf_df.index

ratio_helm  = ratio_helm.rename(columns={'index':'spine'})
ratio_melted_helm["Ref"] = "Helm et al. 2021"
assert (ratio_helm.index == getz_data.index).all()

ratio_merged = pd.concat((ratio_melted_helm,ratio_melted_getz))
compartments = ["surf","int"]
# breakpoint()
for comp in compartments:
    x = "Ref"
    y = "SE"
    xlab = "Method"
    ylab = "Synaptic enrichment"
    labs = []
    pairs = [("Helm et al. 2021","Getz et al. 2022")]
    # comp = "surf"
    fig_file = os.path.join(op_folder,"SE_method_comparison_{0}{1}".format(comp,w_or_wo_ax_label[ax_label]))
    x_pos = [0,1]
    y_pos = np.stack((ratio_helm[comp].to_numpy(),ratio_getz[comp].to_numpy())).T
    plt_widget.SwarmboxLineplotcombo(data=ratio_merged[ratio_merged["localization"]==comp],
                                 x=x,
                                 y=y,
                                 hue=None,
                                 xlab=xlab,
                                 ylab=ylab,
                                 title="",
                                 labs=labs,
                                 pairs=pairs,
                                 color_pal=my_pal4,
                                 stat_test=stat_test,
                                 xfsize=fsize,
                                 yfsize=1.2 * fsize,
                                 fname=fig_file,
                                 save_it=save_it,
                                 ax_lab=ax_label,
                                 x_pos=x_pos,
                                 line_y_data=y_pos)
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    if ax_label == 1:
        xlabel = "SE using Helm et al., 2021"
        ylabel = "SE using Getz et al., 2022"
        title = "{}{}".format(comp[0],subunit)
    else:
        xlabel = "",
        ylabel = "",
        title = ""
    sns.regplot(x=ratio_helm[comp], y=ratio_getz[comp], ax=ax,color=my_pal2[comp])#.\
        # set(xlabel=xlabel,
        #     ylabel=ylabel,
        #     title=title)
    ax.set_xlabel(xlabel,fontsize=fsize)
    ax.set_ylabel(ylabel,fontsize=fsize)
    corr = pearsonr(ratio_helm[comp],ratio_getz[comp])
    fname = os.path.join(op_folder,"SE_correlation_{0}{1}".format(comp,w_or_wo_ax_label[ax_label]))
    ax.text(.05, .8, 'r={:.2f}, p={:.1g}'.format(corr[0],corr[1]),transform=ax.transAxes,fontsize=fsize)
    plt_widget.SaveFigures(fname)
    plt.show()


# breakpoint()