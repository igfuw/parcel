import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest

from parcel import parcel
from libcloudphxx import common as cm

@pytest.fixture(scope="module")  
def data(request):
    SO2_g_init  = 200e-12 
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6 
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12
    outfreq     = 200
    z_max       = 200.
    outfile     = "test_chem_dsl.nc"

    # running parcel model for open chem system  and only for dissolving chem species into droplets
    parcel(dt = .1, z_max = z_max, outfreq = outfreq,\
            SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
            chem_sys = 'open',   outfile = outfile,\
            chem_dsl = True, chem_dsc = False, chem_rct = False,\
            out_bin = \
            '{"radii": {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [3]},\
              "chem" : {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1,\
                  "moms": ["O3_a", "H2O2_a", "SO2_a", "CO2_a", "NH3_a", "HNO3_a"]}}')

    data = netcdf.netcdf_file(outfile,   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", outfile])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "HNO3", "NH3"])
def test_henry_checker(data, chem, eps=2e-5):
                                   #TODO - why so small?
    """                              
    Checking if dissolving chemical compounds into cloud droplets follows Henrys law
    http://www.henrys-law.org/
    """
    conc_aq = data.variables[chem+"_a"][:]
    conc_g  = data.variables[chem+"_g"][:]
    prs     = data.variables["p"][:]
    
    henry = getattr(cm, "H_"+chem), 
    if chem in ["SO2", "CO2", "NH3"]:
        molar_mass = getattr(cm, "M_"+chem+"_H2O")
    else: 
        molar_mass = getattr(cm, "M_"+chem)

    # average drop volume
    vol = np.squeeze(data.variables["radii_m3"][:]) * 4/3. * math.pi

    # dissolved  = partial prsessure * Henry_const * molar mass * drop volume
    henry_aq     = conc_g[1:] * prs[1:] * henry * molar_mass * vol[1:]
    conc_aq      = conc_aq[1:]

    assert np.isclose(conc_aq, henry_aq, atol=0, rtol=eps).all()
