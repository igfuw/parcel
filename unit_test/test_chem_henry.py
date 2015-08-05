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

@pytest.fixture(scope="module")  
def data(request):
    RH_init = .99999
    T_init  = 300.
    p_init  = 100000.
    r_init  = cm.eps * RH_init * cm.p_vs(T_init) / (p_init - RH_init * cm.p_vs(T_init))

    SO2_g_init  = 200e-12 
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6 
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12
    outfreq     = 50
    z_max       = 50.
    outfile     = "test_chem_dsl_"
    dt          = 0.1
    wait        = 1000

    for chem_sys in ["open", "closed"]:
        # running parcel model for open chem system  and only for dissolving chem species into droplets
        parcel(dt = dt, z_max = z_max, outfreq = outfreq, wait=wait,\
                T_0 = T_init, p_0 = p_init, r_0 = r_init,\
                SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
                CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
                chem_sys = chem_sys,   outfile = outfile + chem_sys +".nc",\
                chem_dsl = True, chem_dsc = False, chem_rct = False,\
                out_bin = \
                '{"radii": {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 3]},\
                  "chem" : {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1,\
                      "moms": ["O3_a", "H2O2_a", "SO2_a", "CO2_a", "NH3_a", "HNO3_a"]}}')

    data = {"open":   netcdf.netcdf_file("test_chem_dsl_open.nc",   "r"),\
            "closed": netcdf.netcdf_file("test_chem_dsl_closed.nc", "r")}

    # removing all netcdf files after all tests                                      
    def removing_files():
        for files in data.keys():
            subprocess.call(["rm", "test_chem_dsl_"+files+".nc"])

    request.addfinalizer(removing_files)
    return data

def henry_teor(chem, p, T, vol, conc_g):
    H   = getattr(cm, "H_"  +chem) 
    dHR = getattr(cm, "dHR_"+chem)
    if chem in ["SO2", "CO2", "NH3"]:
        molar_mass = getattr(cm, "M_"+chem+"_H2O")
    else: 
        molar_mass = getattr(cm, "M_"+chem)

    # correction to Henry const. due to temperature
    henry_T = H * np.exp(-1 * dHR * (1./T - 1./298))

    # dissolved  = partial prsessure * Henry_const * molar mass * drop volume
    henry_exp    = conc_g * p * henry_T * molar_mass * vol
    return henry_exp

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "HNO3", "NH3"])
def test_henry_checker(data, chem, eps = 3e-4):
    """                              
    Checking if dissolving chemical compounds into cloud droplets follows Henrys law
    http://www.henrys-law.org/

    libcloudph++ takes into account the affect of temperature on Henry constant and 
    the effects of mass transfer into droplets

    Due o the latter effect to compare with the teoretical values there is first the "wait" period
    with verical velocity set to zero. During this time the droplets may adjust to equilibrum
    and may be then compared with the teoretical values.
    """
    data_open = data["open"]

    vol    = np.squeeze(data_open.variables["radii_m3"][-1]) * 4/3. * math.pi
    T      = data_open.variables["T"][-1]
    p      = data_open.variables["p"][-1]
    conc_g = data_open.variables[chem+"_g"][-1]

    henry_aq     = henry_teor(chem, p, T, vol, conc_g)
    conc_aq      = data_open.variables[chem+"_a"][-1]

    assert np.isclose(conc_aq, henry_aq, atol=0, rtol=eps), chem + " : " + str((conc_aq - henry_aq)/conc_aq) 

def test_henry_plot(data):
    """
    plot mass of dissolved chem species compared with theory
    """
    for chem_sys in ["open", "closed"]:
        plot_henry(data[chem_sys], chem_sys, output_folder = "plots/outputs/")
