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
                      "moms": ["SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                               "CO2_a",  "HCO3_a", "CO3_a", "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'
 
    # run parcel
    parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["SO2", "CO2", "NH3", "HNO3"])
def test_moles_const_dsl_dsc_rct(data, chem, eps =\
                                             {"SO2": 3e-8, "CO2": 2e-14, "NH3": 3e-12, "HNO3":2e-12}):
                                         #TODO why so big?
     """
     Checking if the total number of moles in closed chemical system 
     with dissocoation and chemical reactions remains constant
     
     ini - number of moles in gas phase and aq phase (both ions and not) at t=0
     end - number of moles in gas phase and aq phase (both ions and not) at t=end

     The test is not done for H2O2 and O3 since they are depleted during oxidation
     """
     # NH3 -> NH4+ + OH-
     if chem == "NH3":

         ini = data.variables[chem+"_g"][0] / cm.M_NH3 +\
               data.variables["chem_"+chem+"_a"][0, :].sum() / cm.M_NH3_H2O +\
               data.variables["chem_NH4_a"][0, :].sum()      / cm.M_NH4

         end = data.variables[chem+"_g"][-1] / cm.M_NH3 +\
               data.variables["chem_"+chem+"_a"][-1, :].sum() / cm.M_NH3_H2O +\
               data.variables["chem_NH4_a"][-1, :].sum()      / cm.M_NH4

     # HNO3 -> H+ + NO3-
     if chem == "HNO3":

         ini = data.variables[chem+"_g"][0] / cm.M_HNO3 +\
               data.variables["chem_"+chem+"_a"][0, :].sum() / cm.M_HNO3 +\
               data.variables["chem_NO3_a"][0, :].sum()      / cm.M_NO3

         end = data.variables[chem+"_g"][-1] / cm.M_HNO3 +\
               data.variables["chem_"+chem+"_a"][-1, :].sum() / cm.M_HNO3 +\
               data.variables["chem_NO3_a"][-1, :].sum()      / cm.M_NO3

     # SO2_g -> SO2_a HSO3- SO3-- HSO4- SO4--
     if chem == "SO2":

         ini = data.variables[chem+"_g"][0] / cm.M_SO2 + \
               data.variables["chem_" + chem + "_a"][0, :].sum() / cm.M_SO2_H2O + \
               data.variables["chem_H"+ chem.replace('2','3')+"_a"][0, :].sum() / cm.M_HSO3 + \
               data.variables["chem_" + chem.replace('2','3')+"_a"][0, :].sum() / cm.M_SO3  + \
               data.variables["chem_S_VI"][0, :].sum() / cm.M_H2SO4

         end = data.variables[chem+"_g"][-1] / cm.M_SO2 + \
               data.variables["chem_" + chem + "_a"][-1, :].sum() / cm.M_SO2_H2O + \
               data.variables["chem_H"+ chem.replace('2','3')+"_a"][-1, :].sum() / cm.M_HSO3 + \
               data.variables["chem_" + chem.replace('2','3')+"_a"][-1, :].sum() / cm.M_SO3  + \
               data.variables["chem_S_VI"][-1, :].sum() / cm.M_H2SO4

# TODO
#         print " "
#         print "init = ", data.variables[chem+"_g"][0] / cm.M_SO2, " + ",\
#                          data.variables["chem_" + chem + "_a"][0, :].sum() / cm.M_SO2_H2O, " + ",\
#                          data.variables["chem_H"+ chem.replace('2','3')+"_a"][0, :].sum() / cm.M_HSO3, " + ",\
#                          data.variables["chem_" + chem.replace('2','3')+"_a"][0, :].sum() / cm.M_SO3, " + ",\
#                          data.variables["chem_S_VI"][0, :].sum() / cm.M_H2SO4
#
#         print "init = ", ini
#
#         print " "
#         print "end = ", data.variables[chem+"_g"][-1] / cm.M_SO2, " + ",\
#                          data.variables["chem_" + chem + "_a"][-1, :].sum() / cm.M_SO2_H2O, " + ",\
#                          data.variables["chem_H"+ chem.replace('2','3')+"_a"][-1, :].sum() / cm.M_HSO3, " + ",\
#                          data.variables["chem_" + chem.replace('2','3')+"_a"][-1, :].sum() / cm.M_SO3, " + ",\
#                          data.variables["chem_S_VI"][-1, :].sum() / cm.M_H2SO4
#
#         print "end = ", end
#

     # CO2_g -> CO2_a HCO3- CO3--
     if chem == "CO2":

         ini = data.variables[chem+"_g"][0] / cm.M_CO2 + \
               data.variables["chem_" + chem + "_a"][0, :].sum() / cm.M_CO2_H2O + \
               data.variables["chem_H"+ chem.replace('2','3')+"_a"][0, :].sum() / cm.M_HCO3 + \
               data.variables["chem_" + chem.replace('2','3')+"_a"][0, :].sum() / cm.M_CO3

         end = data.variables[chem+"_g"][-1] / cm.M_CO2 + \
               data.variables["chem_" + chem + "_a"][-1, :].sum() / cm.M_CO2_H2O + \
               data.variables["chem_H"+ chem.replace('2','3')+"_a"][-1, :].sum() / cm.M_HCO3 + \
               data.variables["chem_" + chem.replace('2','3')+"_a"][-1, :].sum() / cm.M_CO3
     
     # do the checking
     assert np.isclose(end, ini, atol=0, rtol=eps[chem]), chem + " : " + str((ini-end)/ini)

def test_H2SO4(data, eps = 3.5e-5):
    """
    Check if the number of dissociated H2SO4 ions is equal to the total number of H2SO4 moles.
    (The library assumes that all H2SO4 dissociates)    

    In libcloudpxx after chemical reactions there is no dissociation step. Therefore, the number
    of HSO4- and SO4-- ions is not the same as the number of H2SO4 moles.
    TODO? 
    
    """

    # read the data
    n_H2SO4 = data.variables["chem_S_VI"][-1, :].sum() / cm.M_H2SO4
    n_HSO4  = data.variables["chem_HSO4_a"][-1, :].sum() / cm.M_HSO4
    n_SO4   = data.variables["chem_SO4_a"][-1, :].sum() / cm.M_SO4

    # number of moles of H2SO4 ions
    n_S6_ions = n_HSO4 + n_SO4

    # do the checking
    assert np.isclose(n_H2SO4, n_S6_ions, atol=0, rtol=eps), str((n_H2SO4 - n_S6_ions) / n_H2SO4)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title='/test_chem_closed_rct_')

