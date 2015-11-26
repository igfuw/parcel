import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")

from parcel import parcel
from libcloudphxx import common
from chem_conditions import parcel_dict

from scipy.io import netcdf
import numpy as np
import pytest
import os, glob
import subprocess
import math

sd_conc_list = [1024, 2*1024, 4*1024, 8*1024, 16*1024, 32*1204, 64*1024, 128*1024, 256*1024]

@pytest.fixture(scope="module")
def data(request):
    """
    run parcel simulations with different initial sd_conc
    return data with values of sd_conc and initial dry mass of the aerosol

    """
    # ... and modify them for the current test
    parcel_dict['chem_dsl'] = True
    parcel_dict['z_max']    = .05
    parcel_dict['outfreq']  = 1
    parcel_dict['out_bin']  = '{"drad": {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [3]}}'

    # lists to store sd_conc at the total dry mass at t=0 of each test run
    out_sd_conc = []
    out_m3_dry  = []

    for sd_conc in sd_conc_list:
        print "sd_conc = ", sd_conc
 
        parcel_dict['outfile'] = "convergence_test_sd_conc=" + str(sd_conc) + ".nc" 
        parcel_dict['sd_conc'] = sd_conc

        # run parcel simulation
        parcel(**parcel_dict)

        # read data
        f_out     = netcdf.netcdf_file(parcel_dict['outfile'], "r")
        mom3_init = f_out.variables["drad_m3"][0,:]
        rhod_init = f_out.variables["rhod"][0]
        chem_rho  = getattr(f_out, "chem_rho")

        # initial dry mass of aerosol [kg/kg dry air]
        ini = mom3_init.sum()  * 4./3 * math.pi * chem_rho
    
        out_sd_conc.append(sd_conc)
        out_m3_dry.append(ini)

    # initial dry mass of aerosol from paper [kg/kg dry air]
    art     = 2.375 * 1e-9 / rhod_init

    data = {"sd_conc" : out_sd_conc, "dry_mass" : out_m3_dry, "art" : art}

    # removing all netcdf files after all tests
    def removing_files():
         for file in glob.glob("convergence_test_sd_conc*"):
            subprocess.call(["rm", file])
    request.addfinalizer(removing_files)
    return data

@pytest.mark.xfail
def test_timestep_print(data):

    print "sd_conc:         ", data["sd_conc"]
    print "dry mass at t=0: ", data["dry_mass"]
    print "article dry mass at t=0: ", data["art"]

