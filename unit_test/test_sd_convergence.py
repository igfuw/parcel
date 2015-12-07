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
import copy

sd_conc_list = [64, 128, 256, 512, 1024, 2*1024, 4*1024, 8*1024, 16*1024, 32*1204]

@pytest.fixture(scope="module")
def data(request):
    """
    run parcel simulations with different initial sd_conc
    return data with values of sd_conc and initial dry mass of the aerosol

    """
    # copy options from chem_conditions
    p_dict = copy.deepcopy(parcel_dict)

    # ... and modify them for the current test
    p_dict['chem_dsl'] = True
    p_dict['z_max']    = .05
    p_dict['outfreq']  = 1
    p_dict['out_bin']  = '{"drad": {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [3]}}'

    # lists to store sd_conc at the total dry mass at t=0 of each test run
    out_sd_conc = []
    out_m3_dry  = []

    for sd_conc in sd_conc_list:
        print "sd_conc = ", sd_conc
 
        p_dict['outfile'] = "convergence_test_sd_conc=" + str(sd_conc) + ".nc" 
        p_dict['sd_conc'] = sd_conc

        # run parcel simulation
        parcel(**p_dict)

        # read data
        f_out     = netcdf.netcdf_file(p_dict['outfile'], "r")
        mom3_init = f_out.variables["drad_m3"][0,:]
        rhod_init = f_out.variables["rhod"][0]
        chem_rho  = getattr(f_out, "chem_rho")

        # initial dry mass of aerosol [kg/kg dry air]
        ini = mom3_init.sum()  * 4./3 * math.pi * chem_rho
    
        out_sd_conc.append(sd_conc)
        out_m3_dry.append(ini)

    data = {"sd_conc" : out_sd_conc, "dry_mass" : out_m3_dry}

    # removing all netcdf files after all tests
    def removing_files():
         for file in glob.glob("convergence_test_sd_conc*"):
            subprocess.call(["rm", file])
    request.addfinalizer(removing_files)
    return data

def test_timestep_print(data, eps=3e-3):
    """
    Check if the total mass of dry aerosol (sum of 3rm moment of dry radii)
    doesn't change too much with different initail super droplet concentration

    """

    dry_mass = np.array(data["dry_mass"]).reshape(data["dry_mass"].__len__());

    # average dry mass from all sd_conc runs
    av = dry_mass.sum() / dry_mass.shape[0]

    for it in range(dry_mass.shape[0]):
        assert np.isclose(av, dry_mass[it], atol=0, rtol=eps), "difference: " + str((av - dry_mass[it]) / av)

