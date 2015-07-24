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
    #TODO - the same for dissociation
    #TODO - think of a similar test for reactions
    SO2_g_init  =  200e-12
    O3_g_init   =  50e-9
    H2O2_g_init =  500e-12
    outfreq     = 20
    z_max       = 200.
    outfile     = "test_chem_closed.nc"
    dt          = .01

    # run parcel
    parcel(dt = dt, z_max = z_max, outfreq = outfreq,\
           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
           chem_sys = 'closed',   outfile = outfile,\
           chem_dsl = True, chem_dsc = False, chem_rct = False,\
           out_bin = '{"chem":\
                        {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": ["O3_a", "H2O2_a", "SO2_a"]}}')

    data = netcdf.netcdf_file(outfile,   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", outfile])

    request.addfinalizer(removing_files)
    return data

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2"])                   #TODO why so big?
def test_chem_closed_dsl(data, chem, eps = {"SO2": 4e-9, "O3": 4e-11, "H2O2": 2e-4}):
    """
    Checking if the total mass of SO_2, O_3 and H2O2 in the closed chemical system 
    with only dissocoation present remains constant

    """
    molar_mass = getattr(cm, "M_"+chem)

    # convert aqueous phase within droplets to gase phase (mole fraction)
    def aq_to_g(chem_aq, rhod, T, M, p):
        return chem_aq * rhod * cm.R * T / M  / p

    ini = data.variables[chem+"_g"][0] 
    #final gas phase = current gas phase + mass dissolved into droplets
    end = data.variables[chem+"_g"][-1] + \
          aq_to_g(data.variables[chem+"_a"][-1], data.variables["rhod"][-1], data.variables["T"][-1],\
                molar_mass, data.variables["p"][-1])

    assert np.isclose(end, ini, atol=0, rtol=eps[chem])
