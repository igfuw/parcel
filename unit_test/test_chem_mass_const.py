import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import pytest
import subprocess
import math
import copy

from parcel import parcel
from libcloudphxx import common as cm
from chem_conditions import *

@pytest.fixture(scope="module")
def data(request):
    """
    Run parcel simulation and return opened netdcf file
    """
    # copy options from chem_conditions ...
    opts_dict = copy.deepcopy(parcel_dict)
    opts_dict['outfreq'] = opts_dict['z_max'] / opts_dict['w'] / opts_dict['dt'] / 4

    # ... and modify them for the current test
    opts_dict['outfile']  = "test_mass.nc"
    opts_dict['chem_dsl'] = True

    opts_dict['out_bin'] = \
             '{"wradii": {"rght": 1e-4, "left": 1e-10, "drwt": "wet", "lnli": "lin", "nbin": 500, "moms": [0, 3]}, \
               "dradii": {"rght": 1e-4, "left": 1e-10, "drwt": "dry", "lnli": "lin", "nbin": 500, "moms": [0, 3]}}'

    # running parcel model
    parcel(**opts_dict)

    data = netcdf.netcdf_file(opts_dict['outfile'], "r")

    # removing all netcdf files after all tests
    def removing_files():
        subprocess.call(["rm", opts_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

def test_water_const(data, eps = 1e-15):
    """
    Check if the total water is preserved

    ini = water vapor mixing ratio at t = 0    + water for aerosol to reach equilibrum at t = 0
    end = water vapor mixing ratio at t = end  + water in all particles (cloud + aerosol) t = end

    """
    rho_w = cm.rho_w 
    rv    = data.variables["r_v"][:]
    mom3  = data.variables["wradii_m3"][:,:]
                                                    
    ini = mom3[0,:].sum()  * 4./3 * math.pi * rho_w + rv[0]
    end = mom3[-1,:].sum() * 4./3 * math.pi * rho_w + rv[-1]

    assert np.isclose(end, ini, atol=0, rtol=eps), str((ini-end)/ini)

def test_dry_mass_const(data, eps = 1e-20):
    """
    Check if the total dry mass is preserved

    ini = dry particulate mass / kg dry air at t=0
    end = dry particulate mass / kg dry air at t=end

    """
    chem_rho = getattr(data, "chem_rho")
    mom3     = data.variables["dradii_m3"][:,:]
    rhod     = data.variables["rhod"][:]
    rv       = data.variables["r_v"][:]

    ini = mom3[0,:].sum()  * 4./3 * math.pi * chem_rho
    end = mom3[-1,:].sum() * 4./3 * math.pi * chem_rho

    assert np.isclose(end, ini, atol=0, rtol=eps), str((ini-end)/ini)
