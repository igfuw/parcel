import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest

from parcel import parcel
from libcloudphxx import common as cm
from chemical_plot import plot_chem

@pytest.fixture(scope="module")
def data(request):

    #TODO - think of a similar test for reactions
    SO2_g_init  = 200e-12
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12

    z_max       = 200.
    dt          = .1
    w           = 1.
    outfreq     = int(z_max / dt / 100)
    sd_conc     = 128.

    outfile     = "test_chem_closed_dsl.nc"

    # run parcel
    parcel(dt = dt, z_max = z_max, outfreq = outfreq, w = w, \
           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
           chem_sys = 'closed',   outfile = outfile,\
           chem_dsl = True, chem_dsc = False, chem_rct = False,\
           sd_conc = sd_conc,\
           out_bin = '{"plt_rw":   {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                       "plt_rd":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                       "plt_ch":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1,\
                                    "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                            "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                            "CO2_a",  "HCO3_a", "CO3_a",\
                                            "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]\
                      }}')

    data = netcdf.netcdf_file(outfile,   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", outfile])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "NH3"])
def test_is_mass_const_dsl(data, chem, eps = {"SO2": 4e-9, "O3": 4e-11, "H2O2": 2e-4, "CO2": 1e-10, "NH3": 2e-7}):
                                                                                #TODO why so big? - check with bigger sd_conc
    """
    Checking if the total mass of SO_2, O_3 and H2O2 in the closed chemical system 
    with only dissolving chem species into droplets, remains constant

    """
    if chem in ["SO2", "CO2", "NH3"]:
        molar_mass = getattr(cm, "M_"+chem+"_H2O")
    else: molar_mass = getattr(cm, "M_"+chem)

    # convert aqueous phase within droplets to gase phase (mole fraction)
    def aq_to_g(chem_aq, rhod, T, M, p):
        return chem_aq * rhod * cm.R * T / M  / p

    ini = data.variables[chem+"_g"][0] 
    #final gas phase = current gas phase + mass dissolved into droplets
    end = data.variables[chem+"_g"][-1] + \
          aq_to_g(data.variables[chem+"_a"][-1], data.variables["rhod"][-1], data.variables["T"][-1],\
                molar_mass, data.variables["p"][-1])

    assert np.isclose(end, ini, atol=0, rtol=eps[chem]), str((ini-end)/ini)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title="/test_chem_closed_dsl_")

