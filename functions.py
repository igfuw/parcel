import math
import numpy as np
from libcloudphxx import common as cm

def mole_frac_to_mix_ratio(mole_frac_g, p, M, T, rhod):
    """
    convert mole fraction [1] to mixing ratio [kg/kg dry air]
    """
    return mole_frac_g * p * M / cm.R / T / rhod

def mix_ratio_to_mole_frac(mix_r, p, M, T, rhod):
    """
    convert mixing ratio [kg/kg dry air] to mole fraction [1] 
    """
    return mix_r * cm.R * T * rhod / p / M

def rh_to_rv(RH, T, p):
    """
    convert relative humidity [%] to water vapor mixing ratio [kg/kg]
    """
    return cm.eps * RH * cm.p_vs(T) / (p - RH * cm.p_vs(T))

def rhod_calc(T, p, rv):
    """
    calculate dry air density
    """
    th_0 = T * (cm.p_1000 / p)**(cm.R_d / cm.c_pd)
    return cm.rhod(p, th_0, rv)

def henry_teor(chem, p, T, vol, mixr_g, rhod):
    """ 
    calculate theoretical mass of chemical species dissolved into cloud droplets - Henry law 
    (per kg of dry air)
    """
    # correction to Henry constant due to temperature
    H       = getattr(cm, "H_"  +chem)
    dHR     = getattr(cm, "dHR_"+chem)
    henry_T = H * np.exp(dHR * (1./T - 1./298))

    # molar mass of chemical species dissolved in water
    if chem in ["SO2", "CO2", "NH3"]:
        M_aq = getattr(cm, "M_"+chem+"_H2O")
    else:
        M_aq = getattr(cm, "M_"+chem)

    # molar mass of chem species in gas phase
    M_g = getattr(cm, "M_"+chem)

    partial_prs = mix_ratio_to_mole_frac(mixr_g, p, M_g, T, rhod) * p

    # dissolved  = partial prsessure * Henry_const * molar mass * drop volume
    return partial_prs * henry_T * M_aq * vol

def henry_teor_2(chem, p, T, vol, mixr_g, rhod):
    """ 
    calculate theoretical mass of chemical species dissolved into cloud droplets - 
    Henry law without temperature correction
    (per kg of dry air)
    """
    H = getattr(cm, "H_"  +chem)

    # molar mass of chemical species dissolved in water
    if chem in ["SO2", "CO2", "NH3"]:
        M_aq = getattr(cm, "M_"+chem+"_H2O")
    else:
        M_aq = getattr(cm, "M_"+chem)

    # molar mass of chem species in gas phase
    M_g = getattr(cm, "M_"+chem)

    partial_prs = mix_ratio_to_mole_frac(mixr_g, p, M_g, T, rhod) * p

    # dissolved  = partial prsessure * Henry_const * molar mass * drop volume
    return partial_prs * H * M_aq * vol

def dissoc_teor(chem, T):
    """ 
    calculate theoretical dissociation constants (as a function of temperature)

    """
    # correction to dissoc constant due to temperature
    K        = getattr(cm, "K_"  +chem)
    dKR      = getattr(cm, "dKR_"+chem)

    print "K =        ", K 
    print "dKR =      ", dKR
    print "exp(...) = ", np.exp(dKR * (1./T - 1./298)) 
    print "K(T) =     ", K * np.exp(dKR * (1./T - 1./298)) 

    return K * np.exp(dKR * (1./T - 1./298))

def log10_size_of_lnr(n_tot, mean_r, lnr, gstdev):
    """
    log-normal size distribution (defined as a function of log_10(radius))

    """
    return n_tot / math.sqrt(2 * math.pi) / math.log(gstdev, 10)\
           * math.exp(-1. * math.pow((lnr - math.log(mean_r, 10)), 2) / 2. / math.pow(math.log(gstdev,10),2))

