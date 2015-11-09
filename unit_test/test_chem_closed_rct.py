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
from kreidenweis_fig1 import plot_fig1
from functions import *

@pytest.fixture(scope="module")
def data(request):

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
    mean_r = .08e-6 / 2
    gstdev = 2./1
    n_tot  = 566.e6

    # process toggling
    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    # output
    z_max       = 1400.
    dt          = .1
    w           = 0.5
    outfreq     = int(z_max / dt / 100) * 4
    sd_conc     = 1024.
    outfile     = "test_chem_closed_rct.nc"

    # run parcel
    parcel(dt = dt, z_max = z_max, outfreq = outfreq, w = w, \
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g = SO2_g_init,  O3_g = O3_g_init,   H2O2_g = H2O2_g_init,\
           CO2_g = CO2_g_init, NH3_g = NH3_g_init, HNO3_g = HNO3_g_init,\
           chem_sys = 'closed',   outfile = outfile,\
           mean_r = mean_r, gstdev = gstdev, n_tot = n_tot, sd_conc=sd_conc, \
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, \
           chem_spn=chem_spn, \
           out_bin = '{\
                  "chem"  : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500,\
                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                      "CO2_a",  "HCO3_a", "CO3_a",\
                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]},\
                  "radii" : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500, "moms": [3]},\
                  "plt_rw": {"rght": 1,    "left": 0,    "drwt": "wet", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_rd": {"rght": 1,    "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_ch": {"rght": 1,    "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,\
                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                      "CO2_a",  "HCO3_a", "CO3_a",\
                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'
    )

    data = netcdf.netcdf_file(outfile,   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", outfile])

    #request.addfinalizer(removing_files)
    return data

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title='/test_chem_closed_rct_')

def test_chem_fig1(data):
    """
    Fig 1 from Kreidenweis et al 2003
    """
    plot_fig1(data, output_folder="plots/outputs", output_title='/Kreidenweis_fig1')

