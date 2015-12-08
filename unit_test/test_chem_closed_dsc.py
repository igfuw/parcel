import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest
import copy

from parcel import parcel
from libcloudphxx import common as cm
from chemical_plot import plot_chem
from chem_conditions import parcel_dict
from functions import dissoc_teor

@pytest.fixture(scope="module")
def data(request):

    # copy options from chem_conditions
    p_dict = copy.deepcopy(parcel_dict)

    # modify options from chem_conditions
    p_dict['sd_conc']  = 1
    p_dict['outfile']  = "test_chem_closed_dsc.nc"
    p_dict['chem_dsl'] = True
    p_dict['chem_dsc'] = True

    p_dict['out_bin']  = p_dict['out_bin'][:-1] + \
        ', "radii" : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 1, "moms": [3]}, \
           "chem"  : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 1,\
                      "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                               "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                               "CO2_a",  "HCO3_a", "CO3_a",\
                               "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'

    # run parcel
    parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    #request.addfinalizer(removing_files)
    return data

def test_is_electroneutral(data, eps = 2e-7):
    """
    Check if after dissociation the electrical charge of cloud droplets is 0

    """
    # read the data
    m_H    = data.variables["chem_H"][-1, :]
    m_OH   = data.variables["chem_OH"][-1, :]
    m_HSO3 = data.variables["chem_HSO3_a"][-1, :]
    m_SO3  = data.variables["chem_SO3_a"][-1, :]
    m_HCO3 = data.variables["chem_HCO3_a"][-1, :]
    m_CO3  = data.variables["chem_CO3_a"][-1, :]
    m_NO3  = data.variables["chem_NO3_a"][-1, :]
    m_NH4  = data.variables["chem_NH4_a"][-1, :]
    m_HSO4 = data.variables["chem_HSO4_a"][-1, :]
    m_SO4  = data.variables["chem_SO4_a"][-1, :]

    # check if positive ions == negative ions 
    # positive ions
    p1 = m_H   / cm.M_H
    p2 = m_NH4 / cm.M_NH4
    n_pos = p1 +p2
    # negative ions
    n1 = m_OH   / cm.M_OH
    n2 = m_HSO3 / cm.M_HSO3 + 2 * m_SO3 / cm.M_SO3
    n3 = m_HCO3 / cm.M_HCO3 + 2 * m_CO3 / cm.M_CO3
    n4 = m_NO3  / cm.M_NO3 
    n5 = m_HSO4 / cm.M_HSO4 + 2 * m_SO4 / cm.M_SO4
    n_neg = n1 + n2 + n3 + n4 + n5

    assert np.isclose(n_neg, n_pos, atol=0, rtol=eps),\
          "electoneutral assert error " + str((n_pos-n_neg)/n_pos) +\
           "\n n_pos = " + str(n_pos) + " and n_neg = " + str(n_neg) + \
           "\n with positive ions: " + \
           "\n     H   " + str(p1) + \
           "\n     NH4 "+ str(p2) + \
           "\n with negative ions: " + \
           "\n     OH  " + str(n1) + \
           "\n     S4  " + str(n2) + \
           "\n     CO2 " + str(n3) + \
           "\n     NO3 " + str(n4) + \
           "\n     S6  " + str(n5)

         
@pytest.mark.parametrize("ion", ["H2O", "SO2", "HSO3", "CO2", "HCO3", "NH3", "HNO3"])
def test_dissoc_constants(data, ion, eps =\
                                {"H2O": 5e-16, "SO2":  8e-6, "HSO3": 6e-6,\
                                 "CO2": 4e-6, "HCO3": 7e-6, "NH3":  2e-6, "HNO3":4e-5}):

     """
     Check if the mass of chemical compounds agrees with the dissociation constants

     For each droplet, from the current mass of chemical compounds in the droplet
     calculate the dissociation constant and compare it with the theoretical value.

     Has to be done per droplet, hence sd_conc = 1 and n_bins = 1

     The accuracy should accually be ~1e-15 
     (that's how it was before adding temperature dependance of dissociation constants).
     However, because the condensation process alters the temperature, there is no way to know
     in what exact temperature the dissociation process takes place. Using the temperature 
     available in parcel output as an approximation of the emperature within droplets during dissociation
     results in lower accuracy.

     """
     # read the data
     V      = data.variables["radii_m3"][-1, :] * 4./3 * math.pi
     T      = data.variables["T"][-1]
     m_H    = data.variables["chem_H"][-1, :]
     m_OH   = data.variables["chem_OH"][-1, :]
     m_SO2  = data.variables["chem_SO2_a"][-1, :]
     m_HSO3 = data.variables["chem_HSO3_a"][-1, :]
     m_SO3  = data.variables["chem_SO3_a"][-1, :]
     m_CO2  = data.variables["chem_CO2_a"][-1, :]
     m_HCO3 = data.variables["chem_HCO3_a"][-1, :]
     m_CO3  = data.variables["chem_CO3_a"][-1, :]
     m_HNO3 = data.variables["chem_HNO3_a"][-1, :]
     m_NO3  = data.variables["chem_NO3_a"][-1, :]
     m_NH4  = data.variables["chem_NH4_a"][-1, :]
     m_NH3  = data.variables["chem_NH3_a"][-1, :]

     # dissociation of water K_H20 = [H][OH]
     def check_water(m_OH, m_H, vol, eps):

         H2O_dissoc = m_OH / cm.M_OH / vol * m_H / cm.M_H / vol

         assert np.isclose(H2O_dissoc, cm.K_H2O, atol=0, rtol=eps),\
                           ion + " : " + str((cm.K_H2O-H2O_dissoc)/cm.K_H2O)

     # dissociation constants K = [A][B]/[AB]
     def check_ions(m_A, M_A, m_B, M_B, m_AB, M_AB, vol, teor_const, eps):

         dissoc_const  = (m_A / M_A * m_B / M_B) / (m_AB / M_AB) / vol
                
         assert np.isclose(dissoc_const, teor_const, atol=0, rtol=eps),\
                           ion + " : " + str((teor_const-dissoc_const)/teor_const)

     # do the checking
     if   ion == "H2O":  check_water(m_OH, m_H, V, eps[ion])
     elif ion == "SO2":\
       check_ions(m_H,   cm.M_H,   m_HSO3, cm.M_HSO3, m_SO2,  cm.M_SO2_H2O, V, dissoc_teor(ion, T),  eps[ion])
     elif ion == "HSO3":\
       check_ions(m_H,   cm.M_H,   m_SO3,  cm.M_SO3,  m_HSO3, cm.M_HSO3,    V, dissoc_teor(ion, T), eps[ion]) 
     elif ion == "CO2":\
       check_ions(m_H,   cm.M_H,   m_HCO3, cm.M_HCO3, m_CO2,  cm.M_CO2_H2O, V, dissoc_teor(ion, T),  eps[ion])
     elif ion == "HCO3":\
       check_ions(m_H,   cm.M_H,   m_CO3,  cm.M_CO3,  m_HCO3, cm.M_HCO3,    V, dissoc_teor(ion, T), eps[ion])
     elif ion == "NH3":\
       check_ions(m_NH4, cm.M_NH4, m_OH,   cm.M_OH,   m_NH3,  cm.M_NH3_H2O, V, dissoc_teor(ion, T),  eps[ion])
     elif ion == "HNO3":\
       check_ions(m_H,   cm.M_H,   m_NO3,  cm.M_NO3,  m_HNO3, cm.M_HNO3,    V, dissoc_teor(ion, T), eps[ion])
     else: assert False

#TODO - add temperature dependance
def test_S6_dissoc(data, eps_HSO4=5e-16, eps_SO4 = 3e-16):
    """
    Check dissociation of H2SO4

    """
    # read the data 
    V      = data.variables["radii_m3"][-1, :] * 4./3 * math.pi
    T      = data.variables["T"][-1]
    m_H    = data.variables["chem_H"][-1, :]
    m_SO4  = data.variables["chem_SO4_a"][-1, :]
    m_HSO4 = data.variables["chem_HSO4_a"][-1, :]
    m_S6   = data.variables["chem_S_VI"][-1, :]

    # dissociation for HSO4 
    left_HSO4 = m_HSO4 / cm.M_HSO4 / V
    rght_HSO4 = (m_H / cm.M_H  / V * m_S6 / cm.M_H2SO4 / V) / (m_H / cm.M_H / V + cm.K_HSO4)

    assert np.isclose(left_HSO4, rght_HSO4 , atol=0, rtol=eps_HSO4),\
                       " HSO4 dissoc error: "+ str((left_HSO4 - rght_HSO4)/left_HSO4)

    # dissociation for SO4
    left_SO4 = m_SO4 / cm.M_SO4 / V
    rght_SO4 = (cm.K_HSO4 * m_S6 / cm.M_H2SO4 / V)  / (m_H / cm.M_H / V + cm.K_HSO4)

    assert np.isclose(left_SO4, rght_SO4 , atol=0, rtol=eps_SO4),\
                       " SO4 dissoc error: " + str((left_SO4 -  rght_SO4)/left_SO4)
 
@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "NH3", "HNO3"])
def test_moles_const_dsl_dsc(data, chem, eps =\
                                {"SO2": 4e-13, "O3":3e-14, "H2O2": 2e-14, "CO2": 2e-14, "NH3": 2e-12, "HNO3":5e-13}):
     """
     Checking if the total number of moles in closed chemical system 
     with dissolving chem species into droplets and dissocoation remains constant

     ini - number of moles in gas phase and aq phase (both ions and non-dissociated) at t=0
     end - number of moles in gas phase and aq phase (both ions and non-dissociated) at t=end

     """
     # O3 and H2O2 (they don't dissociate)
     if  chem in ["O3", "H2O2"] : 
         molar_mass = getattr(cm, "M_"+chem)

         ini = (data.variables[chem+"_g"][0] + data.variables[chem+"_a"][0]) / molar_mass
         end = (data.variables[chem+"_g"][-1] + data.variables[chem+"_a"][-1]) / molar_mass
 
     # NH3 -> NH4+ OH-
     if chem == "NH3":
         ini = data.variables[chem+"_g"][0]                   / cm.M_NH3 +\
               data.variables["chem_"+chem+"_a"][0, :].sum()  / cm.M_NH3_H2O +\
               data.variables["chem_NH4_a"][0, :].sum()       / cm.M_NH4

         end = data.variables[chem+"_g"][-1]                  / cm.M_NH3 +\
               data.variables["chem_"+chem+"_a"][-1, :].sum() / cm.M_NH3_H2O +\
               data.variables["chem_NH4_a"][-1, :].sum()      / cm.M_NH4
 
     # HNO3 -> H+ NO3-
     if chem == "HNO3":
         ini = data.variables[chem+"_g"][0]                   / cm.M_HNO3 +\
               data.variables["chem_"+chem+"_a"][0, :].sum()  / cm.M_HNO3 +\
               data.variables["chem_NO3_a"][0, :].sum()       / cm.M_NO3

         end = data.variables[chem+"_g"][-1]                  / cm.M_HNO3 +\
               data.variables["chem_"+chem+"_a"][-1, :].sum() / cm.M_HNO3 +\
               data.variables["chem_NO3_a"][-1, :].sum()      / cm.M_NO3
 
     # SO2_g -> SO2_a HSO3- SO3-- and CO2 -> CO2_a HCO3- CO3--
     if chem in ["SO2", "CO2"]:
         if chem == "SO2":
             M_gas  = cm.M_SO2
             M_aq   = cm.M_SO2_H2O
             M_ion1 = cm.M_HSO3
             M_ion2 = cm.M_SO3
         elif chem == "CO2":
             M_gas  = cm.M_CO2
             M_aq   = cm.M_CO2_H2O
             M_ion1 = cm.M_HCO3
             M_ion2 = cm.M_CO3
          
         ini = data.variables[chem+"_g"][0]                                     / M_gas + \
               data.variables["chem_"+chem+"_a"][0, :].sum()                    / M_aq + \
               data.variables["chem_H"+chem.replace('2', '3')+"_a"][0, :].sum() / M_ion1 + \
               data.variables["chem_"+chem.replace('2','3')+"_a"][0, :].sum()   / M_ion2

         end = data.variables[chem+"_g"][-1]                                     / M_gas + \
               data.variables["chem_"+chem+"_a"][-1, :].sum()                    / M_aq + \
               data.variables["chem_H"+chem.replace('2', '3')+"_a"][-1, :].sum() / M_ion1 + \
               data.variables["chem_"+chem.replace('2','3')+"_a"][-1, :].sum()   / M_ion2

     # do the checking
     assert np.isclose(end, ini, atol=0, rtol=eps[chem]), chem + " : " + str((ini-end)/ini)

def test_moles_const_S6_dsl_dsc(data, eps=1e-15):
    """
    Check if the number of H2SO4 moles remains constant.
    (There are no chemical reactions, so the mass should stay the same)

    """
    # read the data
    m_S6_ini = data.variables["chem_S_VI"][0, :]
    m_HSO4   = data.variables["chem_HSO4_a"][-1, :]
    m_SO4    = data.variables["chem_SO4_a"][-1, :]
    m_S6     = data.variables["chem_S_VI"][-1, :]

    # convert to moles
    n_S6_ini       = m_S6_ini.sum() / cm.M_H2SO4
    n_S6_end       = m_S6.sum() / cm.M_H2SO4
    n_SO4_HSO4_end = m_SO4.sum() / cm.M_SO4 + m_HSO4.sum() / cm.M_HSO4

    assert np.isclose(n_S6_ini, n_S6_end, atol=0, rtol=eps),       "1: " + str((n_S6_end - n_S6_ini) / n_S6_ini)
    assert np.isclose(n_S6_ini, n_SO4_HSO4_end, atol=0, rtol=eps), "2: " + str((n_SO4_HSO4_end - n_S6_ini) / n_S6_ini)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title='/test_chem_closed_dsc_')

