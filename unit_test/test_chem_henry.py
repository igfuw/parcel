import sys
import os
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
from henry_plot import plot_henry
import functions as fn

@pytest.fixture(scope="module")
def data(request):
    # initial condition
    RH_init = .999999
    T_init  = 300.
    p_init  = 100000.
    r_init  = cm.eps * RH_init * cm.p_vs(T_init) / (p_init - RH_init * cm.p_vs(T_init))

    # calculate rhod for gas init cond
    th_0      = T_init * (cm.p_1000 / p_init)**(cm.R_d / cm.c_pd)
    rhod_init = cm.rhod(p_init, th_0, r_init)

    # init cond for trace gases
    SO2_g_init  = fn.mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
    O3_g_init   = fn.mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
    H2O2_g_init = fn.mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
    CO2_g_init  = fn.mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
    NH3_g_init  = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
    HNO3_g_init = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

    # aerosol size distribution
    mean_r = .08e-6 / 2
    gstdev = 2.
    n_tot  = 566.e6

    # output
    sd_conc     = 1
    outfreq     = 50
    z_max       = 300.
    outfile     = "test_chem_henry.nc"
    dt          = 0.1
    wait        = 100

    parcel(dt = dt, z_max = z_max, outfreq = outfreq, wait=wait,\
          T_0 = T_init, p_0 = p_init, r_0 = r_init,\
          SO2_g = SO2_g_init, O3_g  = O3_g_init,  H2O2_g = H2O2_g_init,\
          CO2_g = CO2_g_init, NH3_g = NH3_g_init, HNO3_g = HNO3_g_init,\
          outfile = outfile,\
          chem_dsl = True, chem_dsc = True, chem_rct = False,\
          aerosol = \
            '{"test": {"kappa": 0.5, "mean_r": [' + str(mean_r) + '], "gstdev": [' + str(gstdev) + '], "n_tot": [' + str(n_tot) + ']}}',\
          sd_conc = sd_conc,\
          out_bin = \
            '{"radii": {"rght": 1.0, "left": 0.0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 3]},\
              "chem" : {"rght": 1.0, "left": 0.0, "drwt": "wet", "lnli": "lin", "nbin": 1,\
               "moms": ["O3_a", "H2O2_a", "SO2_a", "CO2_a", "NH3_a", "HNO3_a", "H"]}}')

    data = netcdf.netcdf_file("test_chem_henry.nc", "r")

    # removing all netcdf files after all tests
    def removing_files():
        subprocess.call(["rm", "test_chem_henry.nc"])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "HNO3", "NH3"])
def test_henry_checker(data, chem, eps = {"SO2": 5e-8, "O3":4e-8, "H2O2": 2e-6, "CO2": 4e-8, "NH3": 4e-7, "HNO3":2e-6}):
    """
    Checking if dissolving chemical compounds into cloud droplets follows Henrys law
    http://www.henrys-law.org/

    libcloudph++ takes into account the effect of temperature and pH on Henry constant
    and the effects of mass transfer into droplets

    Due o the latter effect, to compare with the teoretical values there is first the "wait" period
    with verical velocity set to zero. During this time the droplets may adjust to equilibrum
    and may be then compared with the teoretical values.

    """
    vol    = data.variables["radii_m3"][-1] * 4/3. * math.pi
    conc_H = data.variables["chem_H"][-1, 0] / cm.M_H / vol
    T      = data.variables["T"][-1]
    p      = data.variables["p"][-1]
    mixr_g = data.variables[chem+"_g"][-1]
    rhod   = data.variables["rhod"][-1]

    # dissolved mass according to Henry's law
    henry_aq  = fn.henry_teor(chem, p, T, vol, mixr_g, rhod, conc_H)
    # mass in droplets
    chem_tmp  = data.variables["chem_"+chem+"_a"][-1]

    assert np.isclose(chem_tmp, henry_aq, atol=0, rtol=eps[chem]), chem + " : " + str((chem_tmp - henry_aq)/chem_tmp)

def test_henry_plot(data):
    """
    plot mass of dissolved chem species compared with theory
    """
    plot_henry(data, output_folder = "plots/outputs/")
