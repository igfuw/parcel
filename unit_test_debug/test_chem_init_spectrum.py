# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest
import ast
import copy

from parcel import parcel
from libcloudphxx import common as cm
from chem_conditions import parcel_dict
import functions as fn
from init_spectrum_plot import plot_init_spectrum

@pytest.fixture(scope="module")
def data(request):

    # copy options from chem_conditions
    p_dict = copy.deepcopy(parcel_dict)

    # modify options from chem_conditions
    p_dict['outfile']  = "test_init_spectrum.nc"
    p_dict['chem_dsl'] = True
    p_dict['chem_dsc'] = True

    p_dict['z_max']    = .05
    p_dict['dt']       = .1
    p_dict['w']        = .5
    p_dict['outfreq']  = 1
    p_dict['sd_conc']  = 1024 * 44

    p_dict['out_bin'] = '{"drad": {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [0,3]}}'

    # run parcel
    parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.xfail
def test_init_mass(data, eps = 1e-20):
    """
    Check if the initial mass is the same as in Kreidenweis 2003

    """
    chem_rho = getattr(data, "chem_rho")

    mom3_init = data.variables["drad_m3"][0,:]
    rhod_init = data.variables["rhod"][0]

    # initial dry mass of aerosol [kg/kg dry air]
    ini = mom3_init.sum()  * 4./3 * math.pi * chem_rho

    # initial dry mass of aerosol from paper [kg/kg dry air]
    art = 2.375 * 1e-9 / rhod_init

    print("initial mixing ratio:      ", ini, " [kg / kg dry air]")
    print("article init mixing ratio: ", art, " [kg / kg dry air]")

    assert np.isclose(art, ini, atol=0, rtol=eps), str((ini-art)/ini)

def test_init_spectrum(data, eps = 0.06):
    """
    Checking if the initial aerosol size distribution is close to the analytic solution 

    """
    # size distribution parameters from Kreidenweis 2003
    n_tot   = 566e6
    mean_r  = 0.04e-6
    gstdev  = 2

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
    assert out_bin["drad"]["lnli"] == 'log', "this test should be run with logarithmic spacing of bins"
   
    # parcel initial condition
    rd = data.variables["drad_r_dry"][:] # left bin edges

    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10) 

    # initial size distribution from the model
    model = data.variables['drad_m0'][0,:] * data.variables["rhod"][0] / d_log_rd

    # theoretical solution
    theor     = np.empty(rd.shape) 
    for it in range(rd.shape[0]):                      # evaluate at the middle of the log bin                
        theor[it]   = fn.log10_size_of_lnr(n_tot, mean_r, math.log(rd[it], 10) + d_log_rd / 2, gstdev)

    # check for each bin where there are droplets
    # (TODO missing differences where from empty bin vs analityc solution)
    for it in range(theor.shape[0]):
        if model[it] > 0:
            assert (abs(theor[it] -  model[it]) / model[it] < eps, str(rd[it]) + str(theor[it]) + str(model[it]) + str(abs(theor[it] -  model[it]) / model[it]))

def test_plot_init_spectrum(data):
    """
    Plot the initial dry diameter distribution and compare it with the analitycal solution

    """
    plot_init_spectrum(data, outfolder = "plots/outputs/")
