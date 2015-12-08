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
import pprint 

from parcel import parcel
from libcloudphxx import common as cm
from chemical_plot import plot_chem
from chem_conditions import parcel_dict

@pytest.fixture(scope="module")
def data(request):

    # copy options from chem_conditions
    p_dict = copy.deepcopy(parcel_dict)

    # modify options from chem_conditions
    p_dict['outfile']  = "test_chem_closed_dsl.nc"
    p_dict['chem_dsl'] = True

    pprint.pprint(p_dict)

    # run parcel
    parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "NH3", "HNO3"])
def test_moles_const_dsl(data, chem, eps = {"SO2": 3e-14, "O3": 2e-14,  "H2O2": 2e-14,\
                                            "CO2": 3e-14, "NH3": 2e-14, "HNO3": 2e-14}):
    """
    Checking if the total number of moles in closed chemical system 
    with only dissolving chem species into droplets, remains constant

    ini - number of moles in gas and aq phase at t=0
    end - number of moles in gas and aq phase at t=end

    """
    if chem in ["O3", "H2O2", "HNO3"]:
      M_gas = getattr(cm, "M_"+chem)
      M_aq  = M_gas
    elif chem in ["SO2", "CO2", "NH3"]:
      M_gas = getattr(cm, "M_"+chem)
      M_aq  = getattr(cm, "M_"+chem+"_H2O")

    ini = data.variables[chem+"_g"][0]  / M_gas + data.variables[chem+"_a"][0]  / M_aq
    end = data.variables[chem+"_g"][-1] / M_gas + data.variables[chem+"_a"][-1] / M_aq

    assert np.isclose(end, ini, atol=0, rtol=eps[chem]), chem + " " + str((ini-end)/ini)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title="/test_chem_closed_dsl_")

