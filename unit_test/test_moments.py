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
import pprint as pp

from parcel import parcel
from libcloudphxx import common as cm
from autoconversion import plot_spectrum_m0, plot_profiles

import functions as fn

@pytest.fixture(scope="module")
def data(request):

    p_dict = {}
    p_dict['outfile']  = "test_moms.nc"
    p_dict['sd_conc']  = 1000
    p_dict['outfreq']  = 1
    p_dict['dt']       = 1.
    p_dict['w']        = 0.5
    p_dict['z_max']    = 1.
    p_dict['RH_0']     = 0.9999
    p_dict['T_0']      = 273.15 + 10
    p_dict['p_0']      = 90000.
    p_dict['wait']     = 0

    # initial aerosol: 1-mode ammonium sulfate lognormal
    p_dict['aerosol'] = '{"ammonium_sulfate": {"kappa":  0.61, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot":  [100e6]}}'

    # output for aerosol moments
    p_dict['out_bin'] = '{"aradii": {"rght": 1, "left": 1e-40, "drwt": "dry", "lnli": "log", "nbin": 1, "moms": [0,1,3,6]}}'

    # run parcel
    parcel(**p_dict)

    #simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    #request.addfinalizer(removing_files)
    return data

def test_moments(data):
    """
    test if the initial aerosol size distribution moments agree with the analytical solution for the lognormal size distribution
    """
    import ast
    aerosol = ast.literal_eval(data.aerosol)

    mean   = aerosol['ammonium_sulfate']["mean_r"][0]  # mean radius
    gstdev = aerosol['ammonium_sulfate']["gstdev"][0]  # geometric standard deviation
    n_tot  = aerosol['ammonium_sulfate']["n_tot"][0]   # concentration under standard conditions (T=15C, p=1013.25 hPa, rv=0) [m^-3]

    final_N = n_tot / cm.rho_stp  # calculate the total concentration in [kg^-1 of dry air]
    # moments calculated by the library
    m0 = data.variables["aradii_m0"][0]
    m1 = data.variables["aradii_m1"][0]
    m3 = data.variables["aradii_m3"][0]
    m6 = data.variables["aradii_m6"][0]

    # analytical solution for lognormal size distribution
    m1_a = final_N * mean * math.exp(1./2 * math.pow(math.log(gstdev), 2))
    m3_a = final_N * math.pow(mean, 3) * math.exp(9./2 * math.pow(math.log(gstdev), 2))
    m6_a = final_N * math.pow(mean, 6) * math.exp(36./2 * math.pow(math.log(gstdev), 2))

    assert np.isclose(m0, final_N, atol=0, rtol=1e-3), " m0 diff = " + str((m0 - final_N)/final_N)
    assert np.isclose(m1, m1_a,    atol=0, rtol=1e-3), " m1 diff = " + str((m1 - m1_a)/m1_a)
    assert np.isclose(m3, m3_a,    atol=0, rtol=1e-4), " m3 diff = " + str((m3 - m3_a)/m3_a)
    assert np.isclose(m6, m6_a,    atol=0, rtol=1e-2), " m6 diff = " + str((m6 - m6_a)/m6_a)


