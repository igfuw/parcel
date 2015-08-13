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

    # gas phase
    SO2_g_init  = 200e-12
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12

    # water vapour
    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = cm.eps * RH_init * cm.p_vs(T_init) / (p_init - RH_init * cm.p_vs(T_init))

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
    z_max       = 2400.
    dt          = .1
    w           = 1.
    outfreq     = int(z_max / dt / 100)
    sd_conc     = 1024.
    outfile     = "test_chem_closed_rct.nc"

    # run parcel
#    parcel(dt = dt, z_max = z_max, outfreq = outfreq,\
#           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
#           SO2_g_0 = SO2_g_init,  O3_g_0 = O3_g_init,   H2O2_g_0 = H2O2_g_init,\
#           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
#           chem_sys = 'closed',   outfile = outfile,\
#           mean_r = mean_r, gstdev = gstdev, n_tot = n_tot, sd_conc=sd_conc, \
#           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn=chem_spn, \
#           out_bin = '{\
#                  "chem"  : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500,\
#                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
#                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
#                                      "CO2_a",  "HCO3_a", "CO3_a",\
#                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]},\
#                  "radii" : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500, "moms": [3]},\
#                  "plt_rw": {"rght": 1,    "left": 0,    "drwt": "wet", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
#                  "plt_rd": {"rght": 1,    "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
#                  "plt_ch": {"rght": 1,    "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,\
#                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
#                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
#                                      "CO2_a",  "HCO3_a", "CO3_a",\
#                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'
#    )
#
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

