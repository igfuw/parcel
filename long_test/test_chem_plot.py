import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import pytest
import math
import subprocess

from parcel import parcel
from libcloudphxx import common
from chemical_plot import plot_chem

@pytest.fixture(scope="module")
def data(request):
    # initial condition
    SO2_g_init  = 200e-12
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12
    z_max       = 200. # TODO 2000.
    dt          = .06
    w           = 1.
    outfreq     = int(z_max / dt / 100)
    sd_conc     = 2048.
    
    # turn on chemistry
    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    # define output for moments and chemistry
    out_bin_chem = '{"plt_rw":   {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_rd":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_ch":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1,\
                                  "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                           "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                           "CO2_a",  "HCO3_a", "CO3_a",\
                                           "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]\
                    }}'

    out_bin      = '{"plt_rw": {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_rd": {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]}}'

    # running parcel model for open / closed / off chem system
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
            SO2_g = SO2_g_init,  O3_g = O3_g_init,  H2O2_g = H2O2_g_init,\
            CO2_g = CO2_g_init, NH3_g = NH3_g_init, HNO3_g = HNO3_g_init,\
            chem_sys = 'open',   outfile="test_plot_chem_open.nc",\
            sd_conc = sd_conc,\
            chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
            out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
           SO2_g = SO2_g_init,  O3_g = O3_g_init,  H2O2_g = H2O2_g_init,\
           CO2_g = CO2_g_init, NH3_g = NH3_g_init, HNO3_g = HNO3_g_init,\
           chem_sys = 'closed', outfile="test_plot_chem_closed.nc",\
           sd_conc = sd_conc,\
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
           out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, outfile="test_plot_chem_off.nc",\
           SO2_g=0, O3_g=0, H2O2_g=0, out_bin = out_bin, sd_conc = sd_conc)

    # TODO - why do I have to repeat this import here?
    from scipy.io import netcdf

    data = {'open'   : netcdf.netcdf_file("test_plot_chem_open.nc",   "r"),\
            'closed' : netcdf.netcdf_file("test_plot_chem_closed.nc", "r"),\
            'off'    : netcdf.netcdf_file("test_plot_chem_off.nc",    "r")}

    def removing_files():
        for name, netcdf in data.iteritems():
            netcdf.close()
            subprocess.call(["rm", "test_plot_chem_" + name + ".nc"])

    request.addfinalizer(removing_files)

    return data

def test_chem_plot(data):
    """
    checking if plot function works correctly
    returns quicklook for chemistry
    """
    plot_chem(data, output_folder="plots/outputs", output_title="/test_plot_chem_")
