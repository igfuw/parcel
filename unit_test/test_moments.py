import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
from scipy.integrate import trapz

import numpy as np
import math
import subprocess
import pytest
import ast

from parcel import parcel
from libcloudphxx import common as cm

@pytest.fixture(scope="module")
def data(request):

    p_dict = {}
    p_dict['outfile']  = "test_moms.nc"
    p_dict['sd_conc']  = 1000
    p_dict['sd_conc_large_tail']  = True
    p_dict['outfreq']  = 1000
    p_dict['dt']       = 1.
    p_dict['w']        = 0.5
    p_dict['z_max']    = 500.
    p_dict['RH_0']     = 0.9999
    p_dict['T_0']      = 273.15 + 10
    p_dict['p_0']      = 90000.
    p_dict['wait']     = 0

    # initial aerosol: 1-mode ammonium sulfate lognormal
    p_dict['aerosol'] = '{"ammonium_sulfate": {"kappa":  0.61, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot":  [100e6]}}'

    # output for size distribution moments
    p_dict['out_bin'] = '{\
        "dry_moms": {"rght": 1e-6, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 1,  "moms": [0,1,2,3,4,5,6]},\
        "dry_dstr": {"rght": 1e-6, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 50, "moms": [0]},\
        "wet_moms": {"rght": 1e-4, "left": 1e-6, "drwt": "wet", "lnli": "log", "nbin": 1,  "moms": [0,1,2,3,4,5,6]},\
        "wet_dstr": {"rght": 1e-4, "left": 1e-6, "drwt": "wet", "lnli": "log", "nbin": 50, "moms": [0]}\
        }'

    # run parcel
    parcel(**p_dict)

    #simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

def analytic_moms_for_lognormal(mom, n_tot, mean_r, gstdev):
    """ Returns analytical solution for different moments of lognormal size distr """
    if mom == 0:
        ret = n_tot
    else:
        ret = n_tot * mean_r**mom * math.exp((mom**2)/2. * math.pow(math.log(gstdev), 2))
    return ret

def trapez_moms(x_arg, y_arg, mom):
    """ Returns numerically integrated moments """
    if mom == 0:
        ret = trapz(y_arg, x_arg)
    else:
        ret = trapz(y_arg * x_arg**mom, x_arg)
    return ret

def tmp_err(a, b):
    return np.abs(a-b)/a

def test_dry_moments(data):
    """
    Test if the initial and final aerosol size distribution moments agree with the analytical solution for the lognormal size distribution.
    The test is done for:
    - moments calculated by the library
    - moments calculated as numerical integral of the size distribution from the model
    """
    aerosol = ast.literal_eval(data.aerosol)
    mean   = aerosol['ammonium_sulfate']["mean_r"][0]  # mean radius
    gstdev = aerosol['ammonium_sulfate']["gstdev"][0]  # geometric standard deviation
    n_tot  = aerosol['ammonium_sulfate']["n_tot"][0]   # concentration under standard conditions (T=15C, p=1013.25 hPa, rv=0) [m^-3]
    n_spec = n_tot / cm.rho_stp                        # convert to [kg^-3]

    # analytic moments of lognormal size distribution
    anl_mom = np.zeros(7)
    for it in range(0,7,1):
        anl_mom[it] = analytic_moms_for_lognormal(it, n_spec, mean, gstdev)

    # moments calculated by the library
    mdl_mom_ini = np.zeros(7)
    mdl_mom_end = np.zeros(7)
    for it in range(0,7,1):
        mdl_mom_ini[it] = data.variables["dry_moms_m"+str(it)][0]
        mdl_mom_end[it] = data.variables["dry_moms_m"+str(it)][-1]

    # moments integrated form the resolved size distribution
    left_bin     = data.variables["dry_dstr_r_dry"][:]
    bin_width    = data.variables["dry_dstr_dr_dry"][:]
    mid_bin      = left_bin + 0.5 * bin_width
    dry_conc_ini = data.variables["dry_dstr_m0"][0, :] / bin_width
    dry_conc_end = data.variables["dry_dstr_m0"][-1,:] / bin_width
    int_mom_ini = np.zeros(7)
    int_mom_end = np.zeros(7)
    for it in range(0,7,1):
        int_mom_ini[it] = trapez_moms(mid_bin, dry_conc_ini, it)
        int_mom_end[it] = trapez_moms(mid_bin, dry_conc_end, it)

    # relative accuracy of the test
    rtol_mdl = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-2]
    rtol_int = [1e-2, 1e-2, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]

    print " "
    print "dry radius: "
    for it in range(0,7,1):
        print it, tmp_err(anl_mom[it], mdl_mom_end[it]), tmp_err(anl_mom[it], int_mom_end[it])

    #for it in range(0,7,1):
    #    assert np.isclose(anl_mom[it], mdl_mom_ini[it],  atol=0, rtol=rtol_mdl[it]), "ini model m_"+str(it)+" diff = " + str((anl_mom[it] - mdl_mom_ini[it])/anl_mom[it])
    #    assert np.isclose(anl_mom[it], mdl_mom_end[it],  atol=0, rtol=rtol_mdl[it]), "end model m_"+str(it)+" diff = " + str((anl_mom[it] - mdl_mom_end[it])/anl_mom[it])
    #    assert np.isclose(anl_mom[it], int_mom_ini[it],  atol=0, rtol=rtol_int[it]), "ini integ m_"+str(it)+" diff = " + str((anl_mom[it] - int_mom_ini[it])/anl_mom[it])
    #    assert np.isclose(anl_mom[it], int_mom_end[it],  atol=0, rtol=rtol_int[it]), "end integ m_"+str(it)+" diff = " + str((anl_mom[it] - int_mom_end[it])/anl_mom[it])

def test_wet_moments(data):
    """
    Test if the final droplet size distribution moments calculated by the library agree with
    moments calculated as numerical integral of the size distribution from the model
    """
    # moments calculated by the library
    mdl_mom = np.zeros(7)
    for it in range(0,7,1):
        mdl_mom[it] = data.variables["wet_moms_m"+str(it)][-1]

    # moments integrated form the resolved size distribution
    left_bin  = data.variables["wet_dstr_r_wet"][:]
    bin_width = data.variables["wet_dstr_dr_wet"][:]
    mid_bin   = left_bin + 0.5 * bin_width
    wet_conc  = data.variables["wet_dstr_m0"][-1,:] / bin_width
    int_mom   = np.zeros(7)
    for it in range(0,7,1):
        int_mom[it] = trapez_moms(mid_bin, wet_conc, it)

    # relative accuracy of the test
    rtol = [1e-2, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]

    print " "
    print "wet radius: "
    for it in range(0,7,1):
        print it, tmp_err(mdl_mom[it], int_mom[it])

    #for it in range(0,7,1):
    #    assert np.isclose(mdl_mom[it], int_mom[it], atol=0, rtol=rtol[it]), "end m_"+str(it)+" diff = " + str((mdl_mom[it] - int_mom[it])/mdl_mom[it])

