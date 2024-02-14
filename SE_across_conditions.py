import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib
matplotlib.use("Qt5Agg")
import seaborn as sns
from statannotations.Annotator import Annotator
from SNSPlottingWidget import SNSPlottingWidget
from Utility import COLORS_dict,GluA2StoD, spine_cell_folder,\
    spine_exclude_cells,spine_num_cells, conditions

p_type = "SUM"
stat_test = "Kruskal"

# cell_folder = {"control":"/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/GluA2/Control/stretches/retake",
#               "TTX":"/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/GluA2/TTX/New_data/stretches/TTX",
#               "Untreated":"/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/GluA2/TTX/New_data/stretches/PBS"}
# num_cells = {"TTX":13,"Untreated":10,"control":12}
# exclude_cells = {"TTX":[4,5,6,10,11,12],"Untreated":[1,2],"control":[2,3,4,10]}
All_data = pd.DataFrame()
All_surf_data = pd.DataFrame()
All_int_data = pd.DataFrame()
my_pal1= {"spine": COLORS_dict["spine"],"shaft": COLORS_dict["shaft"]}
my_pal3 = {"internal": COLORS_dict["spine_i"],"surface": COLORS_dict["spine_s"]}
my_pal2 = {"TTX": COLORS_dict["spine"],"control": COLORS_dict["spine_s"],"Untreated":COLORS_dict["shaft"]}

for xdc, condition in enumerate(conditions):
    G2std = GluA2StoD(spine_cell_folder[condition],p_type)
    # excluded_cells =
    int_df,surf_df,gfp_df,\
    int_dend_data,surf_dend_data,\
    gfp_dend_data,int_spine_line_data,\
    surf_spine_line_data,gfp_spine_line_data = G2std.LoadData(spine_num_cells[condition],spine_exclude_cells[condition])

    # gfp_df["dataset"] = 'GFP'
    int_df["condition"] = condition
    surf_df["condition"] = condition
    All_surf_data = pd.concat([All_surf_data,surf_df])
    All_int_data = pd.concat([All_int_data,int_df])
    int_df["dataset"] = 'internal'
    surf_df["dataset"] = 'surface'
    All_data = pd.concat([All_data,int_df])
    All_data = pd.concat([All_data, surf_df])

# All_ratios = pd.DataFrame()
# All_ratios["s_i_ratio"] = All_surf_data[""]
# breakpoint()
num_data=All_data.shape[0]
num_arr = np.arange(0,num_data,1)
All_data.index = num_arr
num_data = All_surf_data.shape[0]
num_arr = np.arange(0,num_data,1)
All_surf_data.index = num_arr
All_int_data.index = num_arr

plt_widget = SNSPlottingWidget()
op_folder = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/GluA2/Figures/Spine_to_dend/"
plt_widget.CreateFolderRecursive(op_folder)

ax_label = 1
save_it = 1
fsize = 22
stat_test = "Kruskal"
ratio_s_d_mean = pd.DataFrame()
stat = "intden"
sp_stat = "sp_{}".format(stat)
sp_bg_stat = "sp_bg_{}".format(stat)
dend_stat = "dend_{}".format(stat)
dend_bg_stat = "dend_bg_{}".format(stat)
ratio_s_d_mean["SE"] = ((All_data[sp_stat]-All_data[sp_bg_stat]))/((All_data[dend_stat] - All_data[dend_bg_stat]))
ratio_s_d_mean["localization"] =  All_data["dataset"]
ratio_s_d_mean["condition"] =  All_data["condition"]
x ="condition"
y="SE"
hue="localization"
palette=my_pal3
labs = ["sGluA2","iGluA2"]
pairs = [
        (("Untreated","surface"),("TTX","surface")),
        (("Untreated","internal"),("TTX","internal"))
        ]
order = ["Untreated","TTX"]
title = ""
xlab = "Condition"
ylab = "Synaptic Enrichment"
hue_order = ["surface","internal"]
fig_file = os.path.join(op_folder,"spine_to_dend_{}_ratio_across_treatment".format(stat))
# breakpoint()
plt_widget.Swarmboxplotcombo(data=ratio_s_d_mean,
                             x=x,
                             y=y,
                             hue=hue,
                             order=order,
                             hue_order=hue_order,
                             xlab= xlab,
                             ylab=ylab,
                             labs=labs,
                             title=title,
                             pairs=pairs,
                             color_pal=my_pal3,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=fsize,
                             fname=fig_file,
                             save_it=save_it)
sp_data = pd.DataFrame()
dend_data = pd.DataFrame()
sp_data[stat] = (All_data[sp_stat]-All_data[sp_bg_stat])
dend_data[stat]  = (All_data[dend_stat]-All_data[dend_bg_stat])
sp_data["localization"] =  All_data["dataset"]
sp_data["condition"] =  All_data["condition"]
dend_data["localization"] =  All_data["dataset"]
dend_data["condition"] =  All_data["condition"]
# breakpoint()
num_bins = 50
print(condition)
fig,ax = plt.subplots(figsize=(12,5),ncols=2,nrows=1)

# ax[1,1].set_title("shaft internal")
cumm = True
stat_to_plot = "probability"
log_scale= True
legend = True
element="poly"
legend = True
alpha = 1
common_norm = False
hue = "condition"
loc = "surface"
fill = False
labs = ["Untreated","TTX"]
ylab = "Cumulative probability"
xlab = "Integrated intesnity"
titles = ["spine {}".format(loc),"shaft {}".format(loc)]
sns.histplot(sp_data[sp_data["localization"] == loc],
             x=stat,
             hue = hue,
             bins=num_bins,
             alpha=alpha,
             stat=stat_to_plot,
             cumulative = cumm,
             legend=True,
             ax=ax[0],
             element=element,
             common_norm=common_norm,
             log_scale=log_scale,
             fill=fill)
sns.histplot(dend_data[dend_data["localization"] == loc],
             x=stat,
             hue=hue,
             hue_order=labs,
             bins=num_bins,
             alpha=alpha,
             stat=stat_to_plot,
             cumulative = cumm,
             legend=legend,
             ax=ax[1],
             element=element,
             log_scale=log_scale,
             common_norm=common_norm,
             fill=fill)
for i,axes in enumerate(ax):
    axes.tick_params(axis='both', which='major', labelsize=18)
    axes.grid(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.legend(labs,frameon=False,fontsize=fsize)
    axes.set_title(titles[i])
if ax_label==1:
    for axes in ax:
        axes.set_xlabel(xlab,fontsize=fsize)
        axes.set_ylabel(ylab,fontsize=fsize)
if save_it == 1:
    plt_widget.SaveFigures(os.path.join(op_folder,"spine_and_dend_surf_{}_hist_across_treatment".format(stat)))
plt.show()

"""
calculating the scaling factor

"""


fig,ax = plt.subplots(figsize=(12,5),ncols=2,nrows=1)

loc = "internal"
titles = ["spine {}".format(loc),"shaft {}".format(loc)]
sns.histplot(sp_data[sp_data["localization"] == loc],
             x=stat,
             hue = hue,
             # hue_order=labs,
             bins=num_bins,
             alpha=alpha,
             stat=stat_to_plot,
             cumulative = cumm,
             legend=True,
             ax=ax[0],
             element=element,
             common_norm=common_norm,
             log_scale=log_scale,
             fill=fill)
sns.histplot(dend_data[dend_data["localization"] == loc],
             x=stat,
             hue=hue,
             # hue_order=labs,
             bins=num_bins,
             alpha=alpha,
             stat=stat_to_plot,
             cumulative = cumm,
             legend=legend,
             ax=ax[1],
             element=element,
             log_scale=log_scale,
             common_norm=common_norm,
             fill=fill)
for i,axes in enumerate(ax):
    axes.tick_params(axis='both', which='major', labelsize=18)
    axes.grid(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    # axes.legend(labs,frameon=False,fontsize=fsize)
    axes.set_title(titles[i])
if ax_label==1:
    for axes in ax:
        axes.set_xlabel(xlab,fontsize=fsize)
        axes.set_ylabel(ylab,fontsize=fsize)
if save_it == 1:
    plt_widget.SaveFigures(os.path.join(op_folder,"spine_and_dend_int_{}_hist_across_treatment".format(stat)))
plt.show()
All_ratios = pd.DataFrame()
stat = "intden"
All_ratios["spine"] = ((All_surf_data['sp_{}'.format(stat)]-All_surf_data['sp_bg_{}'.format(stat)]))/\
                      ((All_int_data['sp_{}'.format(stat)] - All_int_data['sp_bg_{}'.format(stat)]))
All_ratios["shaft"] = ((All_surf_data['dend_{}'.format(stat)]-All_surf_data['dend_bg_{}'.format(stat)]))/\
                      ((All_int_data['dend_{}'.format(stat)] - All_int_data['dend_bg_{}'.format(stat)]))
All_ratios["condition"] = All_surf_data["condition"]
All_ratios["id"] = All_surf_data.index
All_ratios = All_ratios.melt(id_vars=["id","condition"],
                             var_name="localization",
                             value_name="ratio")
# breakpoint()
# fig,ax = plt.subplots()
y="ratio"
hue="localization"
x="condition"
labs = ["spine","shaft"]
pairs = [("TTX","Untreated")]#,("TTX","control"),("control","Untreated")]
order = ["spine","shaft"]

xlab = "Condition"
ylab = r"$\frac{sGlua2}{iGluA2}$"
title = ""
fig_file = os.path.join(op_folder,"surf_to_int_{}_in_spine_across_treatment".format(stat))
plt_widget.Swarmboxplotcombo(data=All_ratios,
                             x=x,
                             y=y,
                             hue=hue,
                             xlab=xlab,
                             ylab=ylab,
                             pairs=pairs,
                             labs=labs,
                             title=title,
                             color_pal=my_pal1,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=1.2*fsize,
                             fname=fig_file,
                             save_it = save_it)

