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

import pprint as pp

from parcel import parcel
from libcloudphxx import common as cm
from chemical_plot import plot_chem
from chem_conditions import parcel_dict
import functions as fn

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
                      "moms": ["O3_a", "H2O2_a", "H", "SO2_a", "S_VI", "CO2_a", "NH3_a", "HNO3_a"]}}'

    pp.pprint(p_dict)

    # run parcel
    parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

def test_is_electroneutral(data, eps = 1e-6):
    """
    Check if after dissociation the electrical charge of cloud droplets is 0
    (check if positive ions == negative ions)

    TODO - why is the accuracy so low?
    """
    # read the data
    m_H  = data.variables["chem_H"][-1, :]
    m_S4 = data.variables["chem_SO2_a"][-1, :]
    m_C4 = data.variables["chem_CO2_a"][-1, :]
    m_N5 = data.variables["chem_HNO3_a"][-1, :]
    m_N3 = data.variables["chem_NH3_a"][-1, :]
    m_S6 = data.variables["chem_S_VI"][-1, :]
    V    = data.variables["radii_m3"][-1, :] * 4./3 * math.pi
    T    = data.variables["T"][-1]

    # helper for concentration of H+
    conc_H = m_H / cm.M_H / V

    # positive ions
    p1 = m_H / cm.M_H
    p2 = fn.diag_n_NH4(m_N3, T, conc_H)

    n_pos = p1 +p2

    # negative ions
    n1 = cm.K_H2O * V / conc_H
    n2 = fn.diag_n_NO3(m_N5, T, conc_H)
    n3 = fn.diag_n_HSO3(m_S4, T, conc_H) + 2 * fn.diag_n_SO3(m_S4, T, conc_H)
    n4 = fn.diag_n_HCO3(m_C4, T, conc_H) + 2 * fn.diag_n_CO3(m_C4, T, conc_H)
    n5 = fn.diag_n_HSO4(m_S6, T, conc_H) + 2 * fn.diag_n_SO4(m_S6, T, conc_H)

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

         
@pytest.mark.parametrize("ion", ["SO2", "HSO3", "CO2", "HCO3", "NH3", "HNO3"])
def test_dissoc_constants(data, ion, eps =\
                                {"SO2": 4e-16, "HSO3": 5e-16,\
                                 "CO2": 4e-16, "HCO3": 5e-16, "NH3":  4e-16, "HNO3":4e-16}):

     """
     Check if the mass of chemical compounds agrees with the dissociation constants

     For each droplet, from the current mass of chemical compounds in the droplet
     calculate the dissociation constant and compare it with the theoretical value.

     Has to be done per droplet, hence sd_conc = 1 and n_bins = 1
     """
     # read the data
     V      = data.variables["radii_m3"][-1, :] * 4./3 * math.pi
     T      = data.variables["T"][-1]

     m_H  = data.variables["chem_H"][-1, :]
     m_S4 = data.variables["chem_SO2_a"][-1, :]
     m_C4 = data.variables["chem_CO2_a"][-1, :]
     m_N5 = data.variables["chem_HNO3_a"][-1, :]
     m_N3 = data.variables["chem_NH3_a"][-1, :]

     # helpers for H+
     conc_H = m_H / cm.M_H / V
     n_H = m_H / cm.M_H

     # dissociation constants K = [A][B]/[AB]
     def check_ions(n_A, n_B, n_AB, vol, teor_const, eps):

         dissoc_const  = n_A * n_B / n_AB / vol
                
         assert np.isclose(dissoc_const, teor_const, atol=0, rtol=eps),\
                           ion + " : " + str((teor_const-dissoc_const)/teor_const)

     # do the checking
     if ion == "SO2":\
       check_ions(n_H, fn.diag_n_HSO3(m_S4, T, conc_H), fn.diag_n_SO2_H2O(m_S4, T, conc_H), V, fn.dissoc_teor(ion, T), eps[ion])
     elif ion == "HSO3":\
       check_ions(n_H, fn.diag_n_SO3(m_S4, T, conc_H),  fn.diag_n_HSO3(m_S4, T, conc_H),    V, fn.dissoc_teor(ion, T), eps[ion]) 
     elif ion == "CO2":\
       check_ions(n_H, fn.diag_n_HCO3(m_C4, T, conc_H), fn.diag_n_CO2_H2O(m_C4, T, conc_H), V, fn.dissoc_teor(ion, T), eps[ion])
     elif ion == "HCO3":\
       check_ions(n_H, fn.diag_n_CO3(m_C4, T, conc_H),  fn.diag_n_HCO3(m_C4, T, conc_H),    V, fn.dissoc_teor(ion, T), eps[ion])
     elif ion == "HNO3":\
       check_ions(n_H, fn.diag_n_NO3(m_N5, T, conc_H),  fn.diag_n_HNO3(m_N5, T, conc_H),    V, fn.dissoc_teor(ion, T), eps[ion])
     elif ion == "NH3":\
       check_ions(fn.diag_n_NH4(m_N3, T, conc_H), cm.K_H2O / conc_H * V, fn.diag_n_NH3_H2O(m_N3, T, conc_H), V, fn.dissoc_teor(ion, T),  eps[ion])
     else: assert False

def test_S6_dissoc(data, eps_HSO4=4e-16, eps_SO4 = 3e-16):
    """
    Check dissociation of H2SO4 

    Done separately because the dissociation formulation assumes no H2SO4 (only ions present)
    and results in different formulation.

    """
    # read the data 
    V      = data.variables["radii_m3"][-1, :] * 4./3 * math.pi
    T      = data.variables["T"][-1]
    m_H    = data.variables["chem_H"][-1, :]
    m_S6   = data.variables["chem_S_VI"][-1, :]

    # helper for H+
    conc_H = m_H / cm.M_H / V

    # dissociation for HSO4 
    left_HSO4 = fn.diag_n_HSO4(m_S6, T, conc_H) / V
    rght_HSO4 = (conc_H * m_S6 / cm.M_H2SO4 / V) / (conc_H + fn.dissoc_teor("HSO4", T))

    assert np.isclose(left_HSO4, rght_HSO4 , atol=0, rtol=eps_HSO4),\
                       " HSO4 dissoc error: "+ str((left_HSO4 - rght_HSO4)/left_HSO4)

    # dissociation for SO4
    left_SO4 = fn.diag_n_SO4(m_S6, T, conc_H) / V
    rght_SO4 = (fn.dissoc_teor("HSO4", T) * m_S6 / cm.M_H2SO4 / V)  / (conc_H + fn.dissoc_teor("HSO4", T))

    assert np.isclose(left_SO4, rght_SO4 , atol=0, rtol=eps_SO4),\
                       " SO4 dissoc error: " + str((left_SO4 -  rght_SO4)/left_SO4)
 
@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "NH3", "HNO3"])
def test_moles_const_dsl_dsc(data, chem, eps =\
                                {"SO2": 2e-14, "O3":3e-14, "H2O2": 5e-15, "CO2": 1.3e-14, "NH3": 3e-15, "HNO3":4e-15}):
     """
     Checking if the total number of moles in closed chemical system 
     with dissolving chem species into droplets and dissocoation remains constant

     ini - number of moles in gas phase and aq phase (both ions and non-dissociated) at t=0
     end - number of moles in gas phase and aq phase (both ions and non-dissociated) at t=end

     """
     # read the data
     n_S4 = data.variables["chem_SO2_a"][:,0]  / cm.M_SO2_H2O
     n_C4 = data.variables["chem_CO2_a"][:,0]  / cm.M_CO2_H2O
     n_N5 = data.variables["chem_HNO3_a"][:,0] / cm.M_HNO3
     n_N3 = data.variables["chem_NH3_a"][:,0]  / cm.M_NH3_H2O

     # do the checking:       
     # O3 and H2O2 (they don't dissociate)
     if  chem in ["O3", "H2O2"] : 
         molar_mass = getattr(cm, "M_"+chem)

         ini = (data.variables[chem+"_g"][0] + data.variables[chem+"_a"][0]) / molar_mass
         end = (data.variables[chem+"_g"][-1] + data.variables[chem+"_a"][-1]) / molar_mass
 
     # NH3 -> NH4+ OH-
     if chem == "NH3":
         ini = data.variables[chem+"_g"][0] / cm.M_NH3 + n_N3[0]
         end = data.variables[chem+"_g"][-1]/ cm.M_NH3 + n_N3[-1]
 
     # HNO3 -> H+ NO3-
     if chem == "HNO3":
         ini = data.variables[chem+"_g"][0] / cm.M_HNO3 + n_N5[0]
         end = data.variables[chem+"_g"][-1]/ cm.M_HNO3 + n_N5[-1]
 
     # SO2_g -> SO2_a HSO3- SO3--
     if chem == "SO2":
         ini = data.variables[chem+"_g"][0] / cm.M_SO2 + n_S4[0]
         end = data.variables[chem+"_g"][-1]/ cm.M_SO2 + n_S4[-1]

     # CO2 -> CO2_a HCO3- CO3--
     if chem == "CO2":
         ini = data.variables[chem+"_g"][0] / cm.M_CO2 + n_C4[0]
         end = data.variables[chem+"_g"][-1]/ cm.M_CO2 + n_C4[-1]

     # do the checking
     assert np.isclose(end, ini, atol=0, rtol=eps[chem]), chem + " : " + str((ini-end)/ini)

def test_moles_const_S6_dsl_dsc(data, eps=1e-15):
    """
    Check if the number of H2SO4 moles remains constant.
    (There are no chemical reactions, so the mass should stay the same)

    """
    # read the data
    m_S6_ini = data.variables["chem_S_VI"][0, :]
    m_S6     = data.variables["chem_S_VI"][-1, :]

    # convert to moles
    n_S6_ini       = m_S6_ini.sum() / cm.M_H2SO4
    n_S6_end       = m_S6.sum() / cm.M_H2SO4

    assert np.isclose(n_S6_ini, n_S6_end, atol=0, rtol=eps),       "1: " + str((n_S6_end - n_S6_ini) / n_S6_ini)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    plot_chem(data, output_folder="plots/outputs", output_title='/test_chem_closed_dsc_')
