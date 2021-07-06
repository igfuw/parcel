import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from math import exp,sqrt, log
from scipy import optimize
import pandas as pd


#Conditon parameters
RH_init = 0.95
T_init  = 288.
p_init  = 100000.

Rd =287;
R_v = 461.5;
c_p = 1004;
eps =0.622;
T_0 = 273.15;
L_v = 2.5e6;
Do = 2.26e-5;
K0 = 2.4e-2
e_s0 = 611.73;
C_gamm = 1.5e-9;
D_cos = 0.9152e-10;
kappak = 0.61;
rho_l = 1000
p_atm = 101325;


##########################################################################################################################

def Supersaturation(RH, press, T):
    p_vs = 6.1078e2 * exp( 17.27 * (T-273.15) / ( T-273.15 + 237.3))
    r_va  = eps * RH * p_vs / (press - RH *p_vs)
    e = press*r_va/(eps+r_va)
    es = e_s0*exp(-((1/(T))-(1/T_0))*(L_v/R_v))
    return (e/es)-1

def pvs(T):
    return 6.1094e2 * exp( 17.625 * (T-273.15) / ( T-273.15 + 243.04))
def lv(T):
    """
    latent heat of evaporation
    """
    return 273.16 + (1850 - 4218) * (T - 273.16)
def lambdaD(T):
    return Do / np.sqrt(2 * R_v * T)
def lambdaK(T, p):
    return (4 / 5) * K0 * T / p / np.sqrt(2 * Rd * T)
def beta(Kn):
    return (1 + Kn) / (1 + 1.71 * Kn + 1.33 * Kn * Kn)
def D(r, T):
    Kn = lambdaD(T) / r  # TODO #57 optional
    return Do * beta(Kn)
def K(r, T, p):
    Kn = lambdaK(T, p) / r
    return K0 * beta(Kn)
def Fd(T, D):
    return 1000 * R_v * T / D / pvs(T)
def Fk(T, K, lv):
    return 1000 * lv / K / T * (lv / R_v / T - 1)
def A(T):
    """
    Koehler curve (expressed in partial pressure)
    """
    return 2 * 0.072 / R_v / T / 1000

def B(kp, rd):
    return kp * rd ** 3

def RH_eq(r, T, kp, rd):
    return 1 + A(T) / r - B(kp, rd) / r ** 3

def dr_dt_MM(r, T, p, RH, kp, rd):
    nom = (RH - RH_eq(r, T, kp, rd))
    den = (
            Fd(T, D(r, T)) +
            Fk(T, K(r, T, p), lv(T))
    )
    return 1 / r * nom / den

def growth2(x, r_d, S, T):
    # return S - (C_gamm/x) + (kappak*pow(r_d,3)/pow(x,3))
    return S - (A(T)/x) + (kappak*pow(r_d,3)/pow(x,3))

def log10_distr(mu, sigma, N, r):

    return (np.exp(-(np.log10(r) - np.log10(mu))**2 / (2 * np.log10(sigma)**2))  / ( np.log10(sigma)* np.sqrt(2 * np.pi)))*N

def log10_distr2(mu, sigma, N, r):

    return (np.exp(-(np.log(r) - np.log(mu))**2 / (2 * np.log(sigma)**2))  / ( np.log(sigma)* np.sqrt(2 * np.pi)))*N*np.log(10)

def r_critical(S):
    licz = 4*C_gamm*C_gamm*C_gamm
    mian = 27*kappak*pow(S-1,2)
    return pow(licz/mian,1/3)
#Properties of log-normal distribution
mu = [0.24e-6, 0.06e-6, 0.54e-6, 0.14e-6]
sigma = [1.58, 1.38, 2.2, 1.2]
N = [90e6, 70e6, 90e6, 80e6]
position = ['left', 'right','left', 'right']

#Setting vectors for calcualtions
radii = np.logspace(-2, 2, 1000) * 1e-6
radii3 = np.logspace(-2, 2, 1000)
radii_wet = np.zeros(len(radii))
radii_2 = np.zeros(len(radii))
drd_drw = np.zeros(len(radii))
x = np.zeros((len(radii)-1, int(len(mu))))
y = np.zeros((len(radii)-1, int(len(mu))))
y2 = np.zeros((len(radii)-1, int(len(mu))))
y3 = np.zeros((len(radii), int(len(mu))))


#Calcualting supersaturation for growth equation
Supe = Supersaturation(RH_init, p_init, T_init)
print(Supe)
def drD_drW(r_w, Sat, T):
    E = exp(C_gamm/r_w)
    A_a = (E - Sat) / (E - (1-kappak)*Sat)
    A_b = pow(A_a,1/3)
    B_a =  C_gamm*kappak*Sat*E/(3*r_w)
    B_b = (E - Sat) * pow((E - (1-kappak)*Sat),2)
    B_c = pow(B_b,-2/3)
    return A_b - B_a * B_c

# Applying growth equation
for i in range(len(radii)):
    # radii_wet[i] = optimize.bisect(growth2, 1e-10, 10, args=(radii[i],Supe-1, T_init))
    radii_2[i] = optimize.bisect(dr_dt_MM, 1e-10, 10, args=(T_init, p_init, RH_init, kappak, radii[i]))
    # drd_drw[i] = drD_drW(np.log10(radii_2[i]), Supe, T_init)
    drd_drw[i] = drD_drW(radii_2[i], 0.95, T_init)* radii_2[i] / radii[i]
print(1+Supe)
# print(drd_drw)
# Plotting
fig1 = plt.figure()
fig1.set_size_inches(18.5, 10.5)
for i in range(len(mu)):
    df = pd.read_csv('lognor_distr_data_'+str(i+1)+'.csv', sep=" ", header=None)
    # Calculating the PDF, pdf(y) = pdf(x) / |dy/dx|
    x[:,i] = log10_distr2(mu[i], sigma[i], N[i], radii[:-1])
    # y[:,i] = log10_distr2(mu[i], sigma[i], N[i], radii[:-1])*(np.diff(np.log10(radii))/np.diff(np.log10(radii_wet)))
    y2[:,i] = log10_distr2(mu[i], sigma[i], N[i], radii[:-1])*(np.diff(np.log(radii))/np.diff(np.log(radii_2)))
    y3[:,i] = log10_distr2(mu[i], sigma[i], N[i], radii)*drd_drw

    # Area under the curve
    pole_w = np.trapz(y[:,i], radii_wet[:-1])
    pole_d = np.trapz(x[:,i], radii[:-1])

    plt.subplot(221+i)
    plt.plot(radii[:-1] * 1e6,  x[:,i] * 1e-6, linestyle=':', label="dry_P")
    plt.plot(df[0] ,  df[1] ,linestyle='dashed', label="dry_H")
    # plt.plot(radii_wet[:-1] * 1e6,  y[:,i] * 1e-6, linestyle='dashed', label="wet_Piotr")
    plt.plot(radii_2[:-1] * 1e6,  y2[:,i] * 1e-6 , linestyle='dashed', label="wet_P")
    plt.plot(radii_2 * 1e6,  y3[:,i] * 1e-6, linestyle=':', label="wet_anlitycznie")
    plt.plot(df[0] ,  df[2] ,linestyle='-.', label="wet_H")

    # plt.plot([r_crit*1e6, r_crit*1e6], [0,100])
    # plt.plot([2, 2], [0, 100])
    plt.xscale('log')
    plt.title("$\mu$= "+str(mu[i]*1e6) + "[μm], $\sigma$ = " + str(sigma[i]) +", N = "+str(N[i]/1e6) +r"[$\frac{\#\cdot 1e6}{m^{-3}}$]",loc=str(position[i]))
    plt.xlabel("particle radius [μm]")
    plt.ylabel(r"$\frac{dN}{dlog_{10}(r)} [cm^{-3}]$")
    plt.xlim((0.01,15))
    # plt.text(0.5, 0.5, str(round(pole_w,4))+'  '+ str(round(pole_d,4)))
    # plt.text(0.4, 1, str((pole_d-pole_w)/pole_d*100)+'%')
    plt.legend()
    plt.grid(True,  linestyle='-.')
    plt.tight_layout(pad=1.9, w_pad=0.5, h_pad=0.2)
fig1.suptitle('Distribution comparison for dry and wet radii for RH ='+str(RH_init)+', T = '+str(T_init)+' [K], P = '+str(p_init/1000)+'[kPa]', fontsize=16)
fig1.savefig('init_spectrum_TEST2.png')
