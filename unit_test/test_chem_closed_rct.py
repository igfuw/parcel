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

@pytest.fixture(scope="module")
def data(request):

    # copy options from chem_conditions
    p_dict = copy.deepcopy(parcel_dict)

    # modify options from chem_conditions
    p_dict['outfile']  = "test_chem_closed_rct.nc"

    p_dict['chem_dsl'] = True
    p_dict['chem_dsc'] = True
    p_dict['chem_rct'] = True

    p_dict['sd_conc'] = 1

    p_dict['out_bin']  = p_dict['out_bin'][:-1] + \
        ', "chem"  : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 1,\
                      "moms": ["SO2_a", "S_VI", "H", "CO2_a", "NH3_a", "HNO3_a", "O3_a", "H2O2_a"]}}'
 
    # run parcel
    parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["CO2", "NH3", "HNO3"])
def test_moles_const_dsl_dsc_rct(data, chem, eps = {"CO2": 2e-14, "NH3": 3e-15, "HNO3":6e-16}):
     """
     Checking if the total number of moles in closed chemical system 
     with dissocoation and chemical reactions remains constant
     
     ini - number of moles in gas phase and aq phase (both ions and not) at t=0
     end - number of moles in gas phase and aq phase (both ions and not) at t=end

     The test is not done for H2O2 and O3 since they are depleted during oxidation
     The test for S4 and S6 moles is done separately below
     """
     # read the data
     n_C4 = data.variables["chem_CO2_a"][:,0]  / cm.M_CO2_H2O
     n_N5 = data.variables["chem_HNO3_a"][:,0] / cm.M_HNO3
     n_N3 = data.variables["chem_NH3_a"][:,0]  / cm.M_NH3_H2O

     # do the checking:       
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

def test_moles_const_S4_S6(data, eps = 1e-14):
     """
     Checking if the total number of moles of S in closed chemical system 
     with dissocoation and chemical reactions remains constant
     
     ini - number of moles in gas phase and aq phase (both ions and not) at t=0
     end - number of moles in gas phase and aq phase (both ions and not) at t=end
     """
     # read the data
     n_S4 = data.variables["chem_SO2_a"][:,0] / cm.M_SO2_H2O
     n_S6 = data.variables["chem_S_VI"][:,0]  / cm.M_H2SO4

     ini_S = data.variables["SO2_g"][0]  / cm.M_SO2 + n_S4[0]  + n_S6[0]
     end_S = data.variables["SO2_g"][-1] / cm.M_SO2 + n_S4[-1] + n_S6[-1]

     # do the checking
     assert np.isclose(end_S, ini_S, atol=0, rtol=eps), "S4 + S6 : " + str((ini_S-end_S)/ini_S)

def test_H2SO4_vs_O3_H2O2(data, eps = 2e-11):
    """
    Check if the increase in S6 moles is equal to the decrease of H2O2 and O3 moles

    """

    ini_O3   = data.variables["O3_g"][0] / cm.M_O3 + \
               data.variables["chem_O3_a"][0, :].sum() / cm.M_O3

    ini_H2O2 = data.variables["H2O2_g"][0] / cm.M_H2O2 + \
               data.variables["chem_H2O2_a"][0, :].sum() / cm.M_H2O2

    ini_S6   = data.variables["chem_S_VI"][0, :].sum() / cm.M_H2SO4


    end_O3   = data.variables["O3_g"][-1] / cm.M_O3 + \
               data.variables["chem_O3_a"][-1, :].sum() / cm.M_O3

    end_H2O2 = data.variables["H2O2_g"][-1] / cm.M_H2O2 + \
               data.variables["chem_H2O2_a"][-1, :].sum() / cm.M_H2O2

    end_S6   = data.variables["chem_S_VI"][-1, :].sum() / cm.M_H2SO4

    # change on O3 and H2O2 moles
    dn_gas = (ini_O3 - end_O3) + (ini_H2O2 - end_H2O2)  

    # change in S6 moles
    dn_s6 = end_S6 - ini_S6

    # do the checking
    assert np.isclose(dn_gas, dn_s6, atol=0, rtol=eps), str((dn_s6 - dn_gas) / dn_gas)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    plot_chem(data, output_folder="plots/outputs", output_title='/test_chem_closed_rct_')

