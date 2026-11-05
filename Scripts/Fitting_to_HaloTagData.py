import json
import os

from GluA2_model_fitting import *
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use("TkAgg")
from SNSPlottingWidget import SNSPlottingWidget
import seaborn as sns
from scipy.stats import spearmanr
sns.set_style(style='white')
from Utility import *

def line_fit(x, m,c):
    return m*x + c

def exp_fit(x,a,b):
    return a*np.exp(-b*x)

def R_seq(y_fit,y_orig):
    ss_res = ((y_orig-y_fit)**2).sum()
    ss_tot = ((y_orig-y_orig.mean())**2).sum()
    # print("in R_Seq =",ss_tot,ss_res)
    return 1 - (ss_res/ss_tot)




in_set = 1  #plot in_set normalized plots ?
save_it = 0 #save the plots or not (=1 for yes)
ax_label = 1 #plot axis labels or not (=1 for yes)
lable_or_no_label = ["","with_label"]
off_set = 0
bin_size = 5
plt_widget = SNSPlottingWidget()
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_folder = os.path.join(repo_root, "Data", "HaloGluA2_analysis")
op_folder = os.path.join(repo_root, "Figures", "HaloGluA2_analysis")
plt_widget.CreateFolderRecursive(op_folder)


# """
## dendritic distribution of the ratio of surface to internal GluA2

to_analyse = [75,100,125,150] #lengths to analyse
mean_ratio = {}
ci_low = {}
ci_high = {}
std_ratio = {}
sem_ratio = {}
mean_surf = {}
std_surf = {}
sem_surf = {}
mean_int = {}
std_int = {}
sem_int = {}
for loi in to_analyse:
    mean_ratio[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_ratio_mean_profile_{}.npy".format(loi)),allow_pickle=True)
    std_ratio[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_ratio_std_profile_{}.npy".format(loi)),allow_pickle=True)
    sem_ratio[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_ratio_sem_profile_{}.npy".format(loi)),allow_pickle=True)
    ci_low[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_ratio_cilow_profile_{}.npy".format(loi)),allow_pickle=True)
    ci_high[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_ratio_cihigh_profile_{}.npy".format(loi)),allow_pickle=True)
    
    mean_surf[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_surf_mean_profile_{}.npy".format(loi)),allow_pickle=True)
    std_surf[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_surf_std_profile_{}.npy".format(loi)),allow_pickle=True)
    sem_surf[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_surf_sem_profile_{}.npy".format(loi)),allow_pickle=True)
    
    mean_int[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_int_mean_profile_{}.npy".format(loi)),allow_pickle=True)
    std_int[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_int_std_profile_{}.npy".format(loi)),allow_pickle=True)
    sem_int[loi] = np.load(os.path.join(data_folder,"Glua2_dendrite_int_sem_profile_{}.npy".format(loi)),allow_pickle=True)

# breakpoint()
fsize = 18
lsize = 14
fig,ax1 = plt.subplots(figsize=(8,6*len(to_analyse)),ncols=1,nrows=len(to_analyse))
if len(to_analyse) ==1:
    ax1 = [ax1]

for ldx,l1 in enumerate(to_analyse):
    # print(ldx,l1)
    x = np.arange(0, l1, bin_size)
    # breakpoint()
    # mean_ratio = mean_ratio[l1].mean(axis=0)[off_set:x.shape[0]+off_set]
    # print("For Length {}, params = ".format(l1),ratio_int_surf[l1].shape,mean_ratio)
    popt, pcov = curve_fit(line_fit, x, mean_ratio[l1])
    popt_e, pcov_e = curve_fit(exp_fit, x, mean_ratio[l1])
    # breakpoint()
    # print("For Length {}, params = ".format(l1),popt_e)
    # print("ratios = ",ratio_int_surf[l1][:,13])
    mean_fit = line_fit(x,*popt)
    mean_fit_exp = exp_fit(x,*popt_e)
    r2 = R_seq(mean_fit,mean_ratio[l1])
    r2_exp = RedChisq(mean_ratio[l1],mean_fit_exp,sem_ratio[l1],len(mean_ratio[l1])-2)
    # print(r2,r2_exp)
    # print("popt = ",popt,"pcov = ",pcov)
    print("popt_e = ",popt_e,"pcov_e = ",pcov_e)
     # ax1[ldx].plot(x,mean_fit,'r-.',label=r"line-fit,$R^2 = %.2f$" % r2)
    ax1[ldx].errorbar(x+bin_size/2,mean_ratio[l1],sem_ratio[l1],color=COLORS_dict["s_2_i_ratio"],label='data',marker='o',
                      markersize=plt_widget.msize,linestyle='None' ,zorder=6)
    ax1[ldx].plot(x+bin_size/2,mean_fit_exp,color=COLORS_dict["s_2_i_ratio"],
                  linestyle='-',label=r"exp-fit" )
    print(r2_exp,1- mean_fit_exp[-1]/mean_fit_exp[0],1-mean_fit[-1]/mean_fit[0])
    ax1[ldx].grid(False)
    # sns.regplot(x,mean_ratio,ax=ax1[ldx],line_kws={"alpha": 0.2})
    pea_ratio = spearmanr(x,mean_ratio[l1])
    # fsize = 18
    loc = "upper right"
    ax1[ldx].set_xlim([-1, l1 + bin_size / 2])
    ax1[ldx].set_ylim([0., 3])
    # ax1[ldx].spines[['left', 'bottom']].set_visible(True)
    # ax1[ldx].yaxis.set_ticks_position('left')
    # ax1[ldx].xaxis.set_ticks_position('bottom')
    ax1[ldx].spines[['right', 'top']].set_visible(False)
    # ax1[ldx].text(.05, .8, 'r={:.2f}, p={:.2g}'.format(pea_ratio[0],pea_ratio[1]),transform=ax1[ldx].transAxes)
    if ax_label == 1:
        ax1[ldx].set_xlabel(r"Dendritic distance ($\mu m$)",size=fsize)
        ax1[ldx].set_ylabel(r"$\frac{sGluA2}{iGluA2}$",size=1.2*fsize)
        # ax1[ldx].set_title("Dendritic distribution of surface to internal fluorescent intesity ratio",size=fsize)
    ax1[ldx].tick_params(axis='both', which='major', labelsize=fsize)
    ax1[ldx].legend(fontsize=lsize,loc=loc,frameon=False)

plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder,"Surface_to_internal_glua2_ratio_exp_fitted"))
plt.show()


fig,ax2 = plt.subplots(figsize=(8,6*len(to_analyse)),ncols=1,nrows=len(to_analyse))
if len(to_analyse) == 1:
    ax2 = [ax2]
# breakpoint()
off_set = 0
for ldx,l1 in enumerate(to_analyse):
    x=  np.arange(0, l1, bin_size)
    ax2[ldx].errorbar(x+bin_size/2,mean_surf[l1],sem_surf[l1],color=COLORS_dict["shaft_s"],label="surface",marker='s',linestyle="-")
    ax2[ldx].errorbar(x+bin_size/2,mean_int[l1],sem_int[l1],color=COLORS_dict["shaft_i"],label="internal",marker='o',linestyle="--")
    # ax2[ldx].errorbar(x+bin_size/2,mean_GFP[l1],sem_GFP[l1],label = "GFP",color="green",marker='8',linestyle="-." )
    # ax2[ldx].set_style(style="white")
    # # ax2[ldx].set_facecolor("white")
    ax2[ldx].spines['right'].set_visible(False)
    ax2[ldx].spines['top'].set_visible(False)
    # ax2[ldx].yaxis.set_ticks_position('left')
    # ax2[ldx].xaxis.set_ticks_position('bottom')
    # ax2[ldx].grid(False)
    if ax_label == 1:
        ax2[ldx].set_ylabel("Normalized fluorescence ",fontsize=fsize)
        ax2[ldx].set_xlabel(r"Dendritic distance ($\mu m$)",fontsize=fsize)#.format(int_glua2_data[l1].shape[0]),fontsize=fsize)
        # ax2[ldx].set_title("Dendritic distribution of normalized fluorescence intensity",fontsize=fsize)
    ax2[ldx].legend(fontsize=lsize,frameon=False,labelcolor='linecolor')
    # ax2[1].legend()
    ax2[ldx].set_xlim([-1,l1+2])
    # ax2[ldx].set_ylim([0,3])
    # ax2[1].set_xlim([-4,100])
plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder,"Glua2_normalized_fluorescent_intensity"))
plt.show()


# """
# Fitting the AMPA GLuA2 model to the normalized fluorescent intensity
# """

ltf = 100 #fitting the 100 micron long data
# mean
x1,ps_ss_fit,pc_ss_fit,ps_chi2,pc_chi2,params_fit,min_func,out_res = FitModel(x,
                                                                              np.stack((mean_surf[ltf],mean_int[ltf])),\
                                                                     np.stack((std_surf[ltf],std_int[ltf])),1.3,np.Inf,ltf)
# breakpoint()
x_grid,ps_ss_dist,pc_ss_dist,pspine_ss_dist = RunSimGluA2(params_fit['dx'].value,
                                                          v_p=10**(params_fit['vp'].value),
                                                          D_s=10**(params_fit['ds'].value),
                                                          D_c=10**(params_fit['dc'].value),
                                                          rat=1.3)

fig,ax = plt.subplots(figsize=(8,6),ncols=1,nrows=1)
ax.tick_params(axis='both', which='major', labelsize=fsize)
plt.plot(x_grid,ps_ss_dist,color=COLORS_dict["shaft_s"],label=r"$P_s$")
plt.plot(x_grid,pc_ss_dist,color=COLORS_dict["shaft_i"],label=r"$P_c$")
plt.plot(x_grid,pspine_ss_dist,color=COLORS_dict["spine_s"],label=r"$P_{spine}$")
plt.title("Surface GluA2")
plt.ylabel("Protein count")
plt.xlabel(r"Distance from soma ($\mu m$)")
plt.legend(fontsize=lsize,labelcolor='linecolor')
plt.show()


fig,ax = plt.subplots(figsize=(8,12),ncols=1,nrows=2)
print("Best fitted D_c =",  10**(params_fit['dc'].value),", D_s =  ",10**(params_fit['ds'].value),", v_p =  ",10**(params_fit['vp'].value))
x = np.arange(0, ltf, bin_size)
ax[0].errorbar(x,mean_surf[ltf],sem_surf[ltf],label = "sGluA2",color=COLORS_dict["shaft_s"],marker='o',
               markersize=plt_widget.msize,linestyle='None' ,zorder=6)
ax[0].plot(x1,ps_ss_fit,c=COLORS_dict["shaft_s"],linestyle="-",label=r"model-fit" )
print(r"surf_fit , $\chi^2_\nu$=%0.2f"%(ps_chi2))
if ax_label==1:
    # ax[0].set_title("Surface GluA2")
    ax[0].set_ylabel("Normalized fluorescence [a.u.]",fontsize=fsize)
    ax[0].set_xlabel(r"Dendritic distance ($\mu m$)",fontsize=fsize)
ax[0].set_ylim([0.2,1.5])
# ax[1].set_ylim([0.2,1.5])
ax[0].set_xlim([-2,ltf+1])
ax[0].tick_params(axis='both', which='major', labelsize=fsize)
ax[0].legend(frameon=False,fontsize=lsize,loc="upper right")
ax[0].spines[['right', 'top']].set_visible(False)
ax[1].errorbar(x,mean_int[ltf],sem_int[ltf],label = "iGluA2",color=COLORS_dict["shaft_i"],marker='o',
               markersize=plt_widget.msize,linestyle='None' ,zorder=6)
ax[1].plot(x1,pc_ss_fit,c=COLORS_dict["shaft_i"],linestyle="-",label=r"model-fit")
ax[1].tick_params(axis='both', which='major', labelsize=fsize)
print(r"int_fit , $\chi^2_\nu$=%0.2f"%(pc_chi2))
if ax_label == 1:
    # ax[1].set_title("Intracellular GluA2");
    ax[1].set_ylabel("Normalized fluorescence [a.u.",fontsize=fsize)
    ax[1].set_xlabel(r"Dendritic distance ($\mu m$)",fontsize=fsize)
ax[1].set_ylim([0.2,1.5])
ax[1].set_xlim([-2,ltf+1])
ax[1].legend(frameon=False,fontsize=lsize,loc="upper right")
ax[1].spines[['right', 'top']].set_visible(False)
fig.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder,"Glua2_normalized_fluorescent_intensity_fitted_len_{}".format(ltf)))
plt.show()
# breakpoint()
# """

