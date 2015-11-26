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
from kreidenweis import plot_fig1
from kreidenweis import plot_fig2
from functions import *
from chem_conditions import *

@pytest.fixture(scope="module")
def data(request):

    # copy options from chem_conditions ...
    opts_dict = copy.deepcopy(parcel_dict)

    # ... and modify them for the current test
    opts_dict['outfile']  = "test_chem_closed_rct.nc"

    opts_dict['chem_dsl'] = True
    opts_dict['chem_dsc'] = True
    opts_dict['chem_rct'] = True
    opts_dict['chem_spn'] = 10

    opts_dict['sd_conc']  = 1024.
    opts_dict['z_max']    = 1400.
    opts_dict['outfreq']  = int(z_max / dt / 100) * 4

    opts_dict['out_bin']  = '{\
                  "chem"  : {"rght": 1e-4, "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 500,\
                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                      "CO2_a",  "HCO3_a", "CO3_a",\
                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]},\
                   "specw": {"rght": 0.0001, "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0, 1, 3]}, \
                   "specd": {"rght": 1e-06,  "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [0, 1, 3]}, \
                  "radii" : {"rght": 1e-4,   "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 500,"moms": [0, 3]},\
                "radiidry": {"rght": 1e-4,   "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 500,"moms": [0, 3]},\
                  "plt_rw": {"rght": 1,      "left": 0,    "drwt": "wet", "lnli": "lin", "nbin": 1,  "moms": [0, 1, 3]},\
                  "plt_rd": {"rght": 1,      "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,  "moms": [0, 1, 3]},\
                  "plt_ch": {"rght": 1,      "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,\
                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                      "CO2_a",  "HCO3_a", "CO3_a",\
                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'
 
    # run parcel
    parcel(**opts_dict)

    #simulation results
    data = netcdf.netcdf_file(opts_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", opts_dict['outfile']])

    request.addfinalizer(removing_files)
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

def test_chem_fig2(data):
    """
    Fig 2 from Kreidenweis et al 2003
    """
    plot_fig2(data, output_folder="plots/outputs", output_title='/Kreidenweis_fig2')



