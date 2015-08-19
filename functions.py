import math
import numpy as np
from libcloudphxx import common as cm

def mole_frac_to_mix_ratio(conc_g, p, M, T, rhod):
    """
    convert mole fraction [1] to mixing ratio [kg/kg dry air]
    """
    return conc_g * p * M / cm.R / T / rhod

def mix_ratio_to_mole_frac(mix_r, p, M, T, rhod):
    """
    convert mixing ratio [kg/kg dry air] to mole fraction [1] 
    """
    return mix_r * cm.R * T * rhod / p / M

def henry_teor(chem, p, T, vol, mixr_g, rhod):
    """ 
    calculate theoretical mass of chemical species dissolved into cloud droplets - Henry law 
    (per kg of dry air)
    """
    H   = getattr(cm, "H_"  +chem)
    dHR = getattr(cm, "dHR_"+chem)
    if chem in ["SO2", "CO2", "NH3"]:
        molar_mass_aq = getattr(cm, "M_"+chem+"_H2O")
    else:
        molar_mass_aq = getattr(cm, "M_"+chem)

    molar_mass_g = getattr(cm, "M_"+chem)

    # correction to Henry const. due to temperature
    henry_T = H * np.exp(-1 * dHR * (1./T - 1./298))

    conc_g = mix_ratio_to_mole_frac(mixr_g, p, molar_mass_g, T, rhod)

    # dissolved  = partial prsessure * Henry_const * molar mass * drop volume
    henry_exp    = conc_g * p * henry_T * molar_mass_aq * vol
    return henry_exp

