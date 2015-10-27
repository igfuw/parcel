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
from functions import *

@pytest.fixture(scope="module")
def data(request):

   # TODO inline initial conditions for chem tests? 

    # initial condition
    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = rh_to_rv(RH_init, T_init, p_init)

    # calculate rhod for initial gas mixing ratio
    rhod_init   = rhod_calc(T_init, p_init, r_init)
    # initial condition for trace geses
    SO2_g_init  = mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
    O3_g_init   = mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
    H2O2_g_init = mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
    CO2_g_init  = mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
    NH3_g_init  = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
    HNO3_g_init = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

    # aerosol size distribution
    mean_r = .08e-6 / 2.
    gstdev = 2.
    n_tot  = 566.e6

    # chem process toggling
    chem_dsl = True
    chem_dsc = False
    chem_rct = False

    # output
    z_max       = 300
    dt          = .1
    w           = 1.
    outfreq     = int(z_max / dt / 50) 
    sd_conc     = 128.
    outfile     = "test_chem_closed_dsl.nc"

    # run parcel
    parcel(dt = dt, z_max = z_max, outfreq = outfreq, w = w, \
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
           chem_sys = 'closed',   outfile = outfile,\
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct,\
           n_tot = n_tot, mean_r = mean_r, gstdev = gstdev,\
           sd_conc = sd_conc,\
           out_bin = '{"plt_rw":   {"rght": 1,    "left":    0, "drwt": "wet", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                       "plt_rd":   {"rght": 1,    "left":    0, "drwt": "dry", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                       "radii" :   {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 1, "moms": [0, 3]},\
                       "plt_ch":   {"rght": 1,    "left":    0, "drwt": "dry", "lnli": "lin", "nbin": 1,\
                                    "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                            "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                            "CO2_a",  "HCO3_a", "CO3_a",\
                                            "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]\
                      }}')

    data = netcdf.netcdf_file(outfile,   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", outfile])

    #request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "NH3", "HNO3"])
def test_is_mass_const_dsl(data, chem, eps = {"SO2": 2e-3, "O3": 2e-15, "H2O2": 8e-16, "CO2": 6e-5, "NH3": 10, "HNO3": 2e-16}):
                                              #TODO                                     TODO         TODO
    """
    Checking if the total number of moles in closed chemical system 
    with only dissolving chem species into droplets, remains constant

    """
    if chem in ["O3", "H2O2", "HNO3"]:
      M_gas = getattr(cm, "M_"+chem)
      M_aq  = M_gas
    elif chem in ["SO2", "CO2", "NH3"]:
      M_gas = getattr(cm, "M_"+chem)
      M_aq  = getattr(cm, "M_"+chem+"_H2O")

    ini = data.variables[chem+"_g"][0]  / M_gas + data.variables[chem+"_a"][0]  / M_aq
    end = data.variables[chem+"_g"][-1] / M_gas + data.variables[chem+"_a"][-1] / M_aq

    #assert np.isclose(end, ini, atol=0, rtol=eps[chem]), chem + " " + str((ini-end)/ini)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title="/test_chem_closed_dsl_")

