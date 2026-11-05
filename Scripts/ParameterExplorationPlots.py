import AMPA_model as AMPM
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
colors = {"surf":'#005f73',"cyto" : '#9b2226',"spine" : '#CA6702'}

class MathTextSciFormatter(mticker.Formatter):
    def __init__(self, fmt="%1.2e"):
        self.fmt = fmt
    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        tup = s.split('e')
        significand = tup[0].rstrip(decimal_point)
        sign = tup[1][0].replace(positive_sign, '')
        exponent = tup[1][1:].lstrip('0')
        if exponent:
            exponent = '10^{%s%s}' % (sign, exponent)
        if significand and exponent:
            s =  r'%s{\times}%s' % (significand, exponent)
        else:
            s =  r'%s%s' % (significand, exponent)
        return "${}$".format(s)

def SaveFigures(filename, ext_list=[".png", ".svg", ".pdf"], dpi=300):
    """
            function to save figures
            required arguments:
                filename
    """

    for ext in ext_list:
        plt.savefig(filename + ext, dpi=dpi)

def ExploreParameter(delta_x,model_L,V_p_list = [],D_s_list=[],D_c_list = [],V_pinit=0,D_sinit=0.1,D_cinit = 1e-3):

    dx = delta_x
    num_rows,num_cols = 3,3
    omega = 60
    fig,ax = plt.subplots(figsize = (30,25),nrows = num_rows,ncols = num_cols)
    x_grid,ps_orig,pc_orig, ps_spine_orig = AMPM.RunSimGluA2(dx,V_pinit,D_cinit,D_sinit,0.4)
    # x_grid = np.arange(0,model_L,dx)
    # f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    # f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    # g = lambda x, pos: "${}$".format(f._formatSciNotation('%1.10e' % x))

    for i in range(num_rows):
        ax[i, 0].plot(x_grid, ps_orig / ps_orig[0],'k--', alpha=1,
                      label=r"orig")
        ax[i, 1].plot(x_grid, pc_orig / pc_orig[0], 'k--', alpha=1 ,
                      label=r"$orig$")
        ax[i, 2].plot(x_grid, ps_spine_orig / omega, 'k--', alpha=1 ,
                  label=r"orig")
        ax[i, 0].set_title(r"$P_s$",fontsize = 14)
        ax[i, 1].set_title(r"$P_c$", fontsize=14)
        ax[i, 2].set_title(r"$P_{spine}$", fontsize=14)
    j = 0
    for i in range(len(V_p_list)):
        x_grid,ps_new,pc_new,ps_spine_new = AMPM.RunSimGluA2(dx,V_p_list[i],D_cinit,D_sinit,0.4)
        # print(mticker.FuncFormatter(V_p_list[i]).)
        ax[j,0].plot(x_grid,ps_new/ps_orig[0],color = colors["surf"],alpha=1/(i*2+1),label=r"$V_p = {:.1E}$".format(V_p_list[i]))
        ax[j,1].plot(x_grid, pc_new/pc_orig[0], color=colors["cyto"], alpha=1 / (i * 2 + 1),label=r"$V_p = {:.1E}$".format(V_p_list[i]))
        ax[j, 2].plot(x_grid, ps_spine_new/omega, color=colors["spine"], alpha=1 / (i * 2+ 1),label=r"$V_p = {:.1E}$".format(V_p_list[i]))
        ax[j,0].set_yscale('log')
        ax[j, 1].set_yscale('log')
        # ax[j, 1].set_yscale('log')
    j = 1
    for i in range(len(D_s_list)):
        x_grid,ps_new,pc_new,ps_spine_new = AMPM.RunSimGluA2(dx,V_pinit,D_cinit,D_s_list[i],0.4)
        ax[j,0].plot(x_grid,ps_new/ps_orig[0],color = colors["surf"],alpha=1/(i*2+1),label=r"$D_s = {:.1E}$".format(D_s_list[i]))
        ax[j,1].plot(x_grid, pc_new/pc_orig[0], color=colors["cyto"], alpha=1 / (i * 2 + 1),label=r"$D_s = {:.1E}$".format(D_s_list[i]))
        ax[j, 2].plot(x_grid, ps_spine_new/omega, color=colors["spine"], alpha=1 / (i * 2 + 1),label=r"$D_s = {:.1E}$".format(D_s_list[i]))
        ax[j, 0].set_yscale('log')
        ax[j, 1].set_yscale('log')
    j = 2
    for i in range(len(D_c_list)):
        x_grid,ps_new,pc_new,ps_spine_new = AMPM.RunSimGluA2(dx,V_pinit,D_c_list[i],D_sinit,0.4)
        ax[j,0].plot(x_grid,ps_new/ps_orig[0],color = colors["surf"],alpha=1/(i*2+1),label=r"$D_c = {:.1E}$".format(D_c_list[i]))
        ax[j,1].plot(x_grid, pc_new/pc_orig[0], color=colors["cyto"], alpha=1 / (i * 2 + 1),label=r"$D_c = {:.1E}$".format(D_c_list[i]))
        ax[j, 2].plot(x_grid, ps_spine_new/omega, color=colors["spine"], alpha=1 / (i * 2 + 1),label=r"$D_c = {:.1E}$".format(D_c_list[i]))
        ax[j, 0].set_yscale('log')
        ax[j, 1].set_yscale('log')

    for i in range(num_rows):
        for k in range(num_cols):
            ax[i,k].legend(fontsize=14)
            ax[i,k].set_xlabel(r"Distance from soma in $\mu m$")
            ax[i, k].set_ylabel(r"Normalized Protein density")
    op_folder = "./"
    fname = "Parameter_exploration"
    plt.tight_layout()
    SaveFigures(os.path.join(op_folder,fname))
    plt.show()


ExploreParameter(delta_x=1,model_L=500,V_p_list = [1e-4,1e-2,1e-1],D_s_list = [0.01,10,100],D_c_list =  [0.01,10,100])

