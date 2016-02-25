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

def henry_teor(chem, p, T, vol, mixr_g, rhod, conc_H):
    """ 
    calculate theoretical mass of chemical species dissolved into cloud droplets - Henry law 
    (per kg of dry air)
    """
    # molar mass of chemical species dissolved in water
    if chem in ["O3", "H2O2"]:
        K1 = 0
        K2 = 0
    elif chem == "SO2":
        K1 = getattr(cm, "K_SO2")
        K2 = getattr(cm, "K_HSO3")
    elif chem == "CO2":
        K1 = getattr(cm, "K_CO2")
        K2 = getattr(cm, "K_HCO3")
    elif chem == "NH3":
        dKR= getattr(cm, "dKR_"+chem)
        K1 = getattr(cm, "K_NH3") * np.exp(dKR * (1./T - 1./298))
        K2 = 0
    elif chem == "HNO3":
        dKR= getattr(cm, "dKR_"+chem)
        K1 = getattr(cm, "K_HNO3") * np.exp(dKR * (1./T - 1./298))
        K2 = 0
    else:
        assert False

    # correction to Henry constant due to temperature and pH
    H       = getattr(cm, "H_"  +chem)
    dHR     = getattr(cm, "dHR_"+chem)
    henry_T = H * np.exp(dHR * (1./T - 1./298)) * (1. + K1/conc_H + K1*K2/conc_H/conc_H)

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

    return K * np.exp(dKR * (1./T - 1./298))

def log10_size_of_lnr(n_tot, mean_r, lnr, gstdev):
    """
    log-normal size distribution (defined as a function of log_10(radius))

    """
    return n_tot / math.sqrt(2 * math.pi) / math.log(gstdev, 10)\
           * math.exp(-1. * math.pow((lnr - math.log(mean_r, 10)), 2) / 2. / math.pow(math.log(gstdev,10),2))

def diag_n_OH(V, conc_H):
    """
    calculate the number of OH moles
    """
    return cm.K_H2O * V / conc_H

def diag_n_NH3_H2O(m_N3, T, conc_H):
    """
    calculate the number of NH3*H2O moles
    """
    return m_N3 / cm.M_NH3_H2O / (1. + dissoc_teor("NH3", T) * conc_H / cm.K_H2O)

def diag_n_NH4(m_N3, T, conc_H):
    """
    calculate the number of NH4+ moles
    """
    return m_N3 / cm.M_NH3_H2O * conc_H * dissoc_teor("NH3", T) / cm.K_H2O / (1. + dissoc_teor("NH3", T) * conc_H / cm.K_H2O)

def diag_n_HNO3(m_N5, T, conc_H):
    """
    calculate the number of HNO3 moles
    """
    return m_N5 / cm.M_HNO3 / (dissoc_teor("HNO3", T) / conc_H + 1)
   
def diag_n_NO3(m_N5, T, conc_H):
    """
    calculate the number of NO3- moles
    """
    return dissoc_teor("HNO3", T) / conc_H * m_N5 / cm.M_HNO3 / (dissoc_teor("HNO3", T) / conc_H + 1)

def diag_n_CO2_H2O(m_C4, T, conc_H):
    """
    calculate the number of CO2*H2O moles
    """
    return m_C4 / cm.M_CO2_H2O \
           / (1 + dissoc_teor("CO2", T) / conc_H + dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H)

def diag_n_HCO3(m_C4, T, conc_H):
    """
    calculate the number of HCO3- moles
    """
    return m_C4 / cm.M_CO2_H2O * dissoc_teor("CO2", T) / conc_H \
           / (1 + dissoc_teor("CO2", T) / conc_H + dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H)

def diag_n_CO3(m_C4, T, conc_H):
    """
    calculate the number of CO3-- moles
    """
    return m_C4 / cm.M_CO2_H2O * dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H \
           / (1 + dissoc_teor("CO2", T) / conc_H + dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H)

def diag_n_SO2_H2O(m_S4, T, conc_H):
    """
    calculate the number of SO2*H2O moles
    """
    return m_S4 / cm.M_SO2_H2O \
           / (1 + dissoc_teor("SO2", T) / conc_H + dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H)

def diag_n_HSO3(m_S4, T, conc_H):
    """
    calculate the number of HSO3- moles
    """
    return m_S4 / cm.M_SO2_H2O * dissoc_teor("SO2", T) / conc_H \
           / (1 + dissoc_teor("SO2", T) / conc_H + dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H)


def diag_n_SO3(m_S4, T, conc_H):
    """
    calculate the number of SO3-- moles
    """
    return m_S4 / cm.M_SO2_H2O * dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H \
           / (1 + dissoc_teor("SO2", T) / conc_H + dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H)

def diag_n_HSO4(m_S6, T, conc_H):
    """
    calculate the number of HSO4- moles
    """
    return m_S6 / cm.M_H2SO4 * conc_H / (conc_H + dissoc_teor("HSO4", T))

def diag_n_SO4(m_S6, T, conc_H):
    """
    calculate the number of SO4-- moles
    """
    return m_S6 / cm.M_H2SO4 * dissoc_teor("HSO4", T) / (conc_H + dissoc_teor("HSO4", T))
