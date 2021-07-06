import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from math import exp,sqrt, log
from scipy import optimize

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

def growth2(x, r_d, S):
    return S - (C_gamm/x) + (kappak*pow(r_d,3)/pow(x,3))

def log10_distr(mu, sigma, N, r):

    return (np.exp(-(np.log10(r) - np.log10(mu))**2 / (2 * np.log10(sigma)**2))  / ( np.log10(sigma)* np.sqrt(2 * np.pi)))*N

mu = [0.24e-6, 0.06e-6, 0.54e-6, 0.14e-6]
sigma = [1.58, 1.38, 2.2, 1.1]
N = [90e6,70e6, 90e6, 80e6]
position = ['left', 'right','left', 'right']


radii = np.logspace(-2, 2, 1000) * 1e-6
radii_wet = np.zeros(len(radii))
radii_2 = np.zeros(len(radii))
drd_drw = np.zeros(len(radii))
Supe = Supersaturation(RH_init, p_init, T_init)

def drD_drW(r_w, S):
    je = 2*C_gamm/3/kappak
    dw = S/kappak
    tr =  C_gamm/kappak
    return (je * pow(r_w,2) - dw * pow(r_w,3)) / pow(tr * pow(r_w,2) - dw * pow(r_w,3),2/3)




for i in range(len(radii)):
    radii_wet[i] = optimize.bisect(growth2, 1e-10, 10, args=(radii[i],Supe))
    radii_2[i] = optimize.bisect(dr_dt_MM, 1e-10, 10, args=(T_init, p_init, RH_init, kappak, radii[i]))
    drd_drw[i] = drD_drW(radii_wet[i], Supe)



x = np.zeros((len(radii)-1, int(len(mu))))
y = np.zeros((len(radii)-1, int(len(mu))))
y2 = np.zeros((len(radii)-1, int(len(mu))))
# y2 = np.zeros((len(radii), int(len(mu))))


fig1 = plt.figure()
fig1.set_size_inches(18.5, 10.5)
for i in range(len(mu)):
    x[:,i] = log10_distr(mu[i], sigma[i], N[i], radii[:-1])
    y[:,i] = log10_distr(mu[i], sigma[i], N[i], radii[:-1])/(np.diff(np.log10(radii_wet))/np.diff(np.log10(radii)))
    y2[:,i] = log10_distr(mu[i], sigma[i], N[i], radii[:-1])/(np.diff(np.log10(radii_wet))/np.diff(np.log10(radii)))
    # y2[:,i] = log10_distr(mu[i], sigma[i], N[i], radii[:])*abs(drd_drw)
    pole_w = np.trapz(y[:,i], radii_wet[:-1])
    pole_d = np.trapz(x[:,i], radii[:-1])
        # pole_w = np.trapz(y[:,i], np.diff(np.log10(radii_wet)))
        # pole_d = np.trapz(x[:,i], np.diff(np.log10(radii)))
    # pole_w = np.trapz(y2[:-1,i], np.diff(radii_wet))
    plt.subplot(221+i)
    # plt.plot(radii[:-1]*1e6, radii_wet[:-1]/radii[:-1])
    # plt.plot(radii_wet[:-1], radii[:-1])
    plt.plot(radii[:-1] * 1e6,  x[:,i] * 1e-6,  label="dry")
    plt.plot(radii_wet[:-1] * 1e6,  y[:,i] * 1e-6,  label="wet_Piotr")
    plt.plot(radii_2[:-1] * 1e6,  y2[:,i] * 1e-6,  label="wet_Sylwester")
    # plt.plot(radii_wet[:] * 1e6,  y2[:,i] * 1e-6,  label="wet_equ")
    plt.xscale('log')
    plt.title("$\mu$= "+str(mu[i]*1e6) + "[m], $\sigma$ = " + str(sigma[i]) +", N = "+str(N[i]/1e6) +r"[$\frac{\#\cdot 1e6}{m^{-3}}$]",loc=str(position[i]))
    plt.xlabel("particle radius [μm]")
    plt.ylabel(r"$\frac{dN}{dlog_{10}(D)} [cm^{-3}]$")
    plt.xlim((0.01,15))
    plt.text(0.5, 0.5, str(round(pole_w,4))+'  '+ str(round(pole_d,4)))
    # plt.text(0.4, 1, str((pole_d-pole_w)/pole_d*100)+'%')
    plt.legend()
    plt.grid(True,  linestyle='-.')
    plt.tight_layout(pad=1.9, w_pad=0.5, h_pad=0.2)
fig1.suptitle('Distribution comparison for dry and wet radii for RH ='+str(RH_init)+', T = '+str(T_init)+' [K], P = '+str(p_init/1000)+'[kPa]', fontsize=16)
fig1.savefig('init_spectrum_TEST2.png')


#lognormal = st.lognorm(s=sigma, scale = np.exp(mu))
#theor = np.empty(radii.shape)
#for it in range(radii.shape[0]):
#    theor[it] = log10_size_of_lnr(N, mu, math.log(radii[it], 10), sigma)
#    x[it] = (np.exp(-(np.log10(radii[it]) - np.log10(mu))**2 / (2 * np.log10(sigma)**2))  / ( np.log10(sigma) * np.sqrt(2 * np.pi)))*N
#    x[it] = st.lognorm.pdf(radii[it],loc=np.log10(mu), scale=np.exp(mu), s=sigma)*N
#print(theor)
#
#
#x = lognormal.pdf(radii)*N
#print(x)
