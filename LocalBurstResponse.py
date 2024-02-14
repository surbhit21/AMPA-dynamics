import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
import numpy as np
import scipy as sc

"""
This file provides the solution for the equation of local translation burst in Eqn 22 of Fonkeu et al. 2019, Neuron Paper
The equation can be interpreted as 
1. Local translation of mRNA such as CNIH2 in response to burst of activity
2. Spread of a local increase in activation of kinases and phosphatases such as CaMK2 and Calcineurine
"""


def Eerf(A, B, C):
    # print("A,B,C = ",A,C)
    return np.exp(2 * A * B) * sc.special.erf(A * C + B / C) + np.exp(-2. * A * B) * sc.special.erf(A * C - B / C)


def Integral_minus1_2(a, b, A, B):
    # print("a,b,A,B = ",a,b, A)
    # a, b, A, B = np.abs(a), np.abs(b), np.abs(A), np.abs(B)
    # breakpoint()
    return 0.5 * np.emath.sqrt(np.pi / A)  * (Eerf(np.emath.sqrt(A) , np.emath.sqrt(B) , np.emath.sqrt(b) ) - Eerf(np.emath.sqrt(A) , np.emath.sqrt(B) , np.emath.sqrt(a) ))


def Integral_1_2(a, b, A, B):
    # print("a,b,A,B = ", a, b, A)
    # a,b,A,B =  np.abs(a),np.abs(b),np.abs(A),np.abs(B)
    # breakpoint()
    return  ((-1. / A) * (np.emath.sqrt(b) * np.exp(-1. * A * b - B / b) - np.emath.sqrt(a) * np.exp(-1. * A * a - B / a))) + \
           ((1. / (4 * A)) * np.emath.sqrt(np.pi / A) * (
                       Eerf(np.emath.sqrt(A) , np.emath.sqrt(B) , np.emath.sqrt(b) ) - Eerf(np.emath.sqrt(A) , np.emath.sqrt(B) , np.emath.sqrt(a) ))) + \
           ((np.emath.sqrt(B * np.pi) / (2 * A)) * (
                       Eerf(np.emath.sqrt(A) , np.emath.sqrt(B), np.emath.sqrt(1 / b)) - Eerf(np.emath.sqrt(A) , np.emath.sqrt(B) , np.sqrt(1 / b) )))
    # return soln

def LocalBurst(t_p, deg_c, x_c, t_c, L, dx):
    Dp, Vp = t_p  # t_p describes the transport of species
    Kp = deg_c  # either degradation rate of protein or deactivation rate of kinases
    beta_b, x_0, sigma = x_c  # controls the amplitude, location and spatial spread
    t_0, a1, a2, b1, b2, T_max = t_c  # controls the on time, amplitude and time constant of burst
    x_grid = np.arange(0, L, dx)
    dt = 0.1
    t_grid = np.arange(0, 40*60, dt)
    P_X_T = np.empty((t_grid.shape[0], x_grid.shape[0]),dtype=np.float64)
    # P_X_T[0] = np.ones(x_grid.shape[0])
    print(P_X_T.shape)
    # some auxiliary variables
    # breakpoint()
    eta_p = Kp + Vp ** 2 / (4. * Dp)
    mu_p = Vp / (2. * Dp)
    a = sigma
    A1 = (eta_p - a2) / (4. * Dp)
    A2 = (eta_p - b2) / (4. * Dp)
    B = (sigma * mu_p + 2 * (x_grid - x_0)) ** 2 / 4.
    print(a,A1,A2,B)
    for idx,i in enumerate(t_grid[int(t_0/dt):]):

        b = sigma + 4. * Dp * (i - t_0)

        I_m12_A1_B = Integral_minus1_2(a, b, A1, B)
        I_12_A1_B = Integral_1_2(a, b, A1, B)
        I_m12_A2_B = Integral_minus1_2(a, b, A2, B)
        I_12_A2_B = Integral_1_2(a, b, A2, B)
        e1 = beta_b * mu_p * mu_p / 2. + sigma * Kp / (4. * Dp) + mu_p * (x_grid - x_0)
        # print(e1.shape,P_X_T[i].shape)
        # breakpoint()
        P_X_T[idx] = ((beta_b * np.sqrt(sigma) * np.exp(e1)) / (16 * T_max * Dp *Dp)) \
                   * ((a1 * np.exp(-a2 * (i - t_0 + sigma / (4 * Dp))) * (
                    I_m12_A1_B * (4 * Dp*(i - t_0) + sigma) - I_12_A1_B)) \
                      + (b1 * np.exp(-b2 * (i - t_0 + sigma / (4 * Dp))) * (
                            I_12_A2_B * (4 * Dp*(i - t_0) + sigma) - I_m12_A2_B)))

    return x_grid,t_grid,P_X_T

def Phi(t,t0):
    if t >= t0:
        return 1
    return 0
def get_temporal_profile(t_c,t_grid):
    t0, a1, a2, b1, b2, T_max = t_c
    t_profile = np.zeros(t_grid.shape)
    for tdx,t in enumerate(t_grid):
        t_profile[tdx] = Phi(t,t0)*(t-t0)*(a1*np.exp(-a2*(t-t0)) - b1*np.exp(-b2*(t-t0)))/T_max
    return t_profile

def get_spatial_profile(x_c,x_grid):
    beta_b, x_0, sigma = x_c
    x_prof = beta_b*np.exp(-((x_grid-x_0)**2)/sigma)
    return x_prof

def get_burst_profile(x_c,x_grid,t_c,t_grid):
    t_profile = get_temporal_profile(t_c,t_grid)
    x_profile = get_spatial_profile(x_c,x_grid)
    breakpoint()
    return np.outer(t_profile,x_profile.T)
tp = [0.1,1e-5]
LifTP = 5.6*24*60*60;
kp = np.log(2) / LifTP;
L = 500
dx = 0.24
betab = 0.21
location = 250
xc = [betab,location,dx]
tc = [0,1,1/(12),0.1,1/(12),4]

XG,TG,PXT = LocalBurst(tp,kp,xc,tc,L,dx)
t_prof = get_temporal_profile(tc,TG)
x_prof = get_spatial_profile(xc,XG)
bp_profile = get_burst_profile(xc,XG,tc,TG)
fig,ax = plt.subplots()
TG /= 60
# XG /= dx
im = ax.contourf(XG,TG,bp_profile,origin="upper")
# divider2 = make_axes_locatable(axa)
# cax2 = divider2.append_axes(pos, size='5%', pad=pad)
plt.colorbar(im, ax=ax);
# plt.show()
# breakpoint()
# TG /= 60
fig,ax = plt.subplots()
plt.plot(TG,t_prof)
# plt.plot(x_prof,XG/dx)
# plt.show()

fig,ax = plt.subplots()
im = ax.contourf(XG,TG,PXT,origin="upper")
# divider2 = make_axes_locatable(axa)
# cax2 = divider2.append_axes(pos, size='5%', pad=pad)
plt.colorbar(im, ax=ax);
plt.show()
# breakpoint()
# x= np.arange(0,10,1)
# y = np.exp(x)
# plt.plot(x,y)
# plt.show()
# P_X_T

