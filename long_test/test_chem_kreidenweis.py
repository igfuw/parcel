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
    p_dict = copy.deepcopy(parcel_dict)

    # ... and modify them for the current test
    p_dict['outfile']  = "test_chem_kreidenweis.nc"

    p_dict['chem_dsl'] = True
    p_dict['chem_dsc'] = True
    p_dict['chem_rct'] = True
    p_dict['chem_spn'] = 10

    p_dict['sd_conc']  = 256 
    p_dict['z_max']    = 1400.
    p_dict['outfreq']  = int(z_max / dt / 100) * 4

    p_dict['out_bin']  = '{\
                  "chem"  : {"rght": 1e-4,  "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 50, "moms": ["H"]},\
                  "radii" : {"rght": 1e-4,  "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 50, "moms": [3]},\
                   "specd": {"rght": 1e-06, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 50,  "moms": [0, 1, 3]},\
                  "plt_rw": {"rght": 1,     "left": 0,     "drwt": "wet", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_rd": {"rght": 1,     "left": 0,     "drwt": "dry", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_ch": {"rght": 1,     "left": 0,     "drwt": "dry", "lnli": "lin", "nbin": 1,\
                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                      "CO2_a",  "HCO3_a", "CO3_a",\
                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'

    # run parcel
    parcel(**p_dict)

    #simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    #request.addfinalizer(removing_files)
    return data

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title='/test_chem_kreidenwies_')

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
