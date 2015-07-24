import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf
import numpy as np
import math
import subprocess

from libcloudphxx import common as cm
from parcel import parcel

def test_chem_closed():
    """
    Checking if the total mass of SO_2, O_3 and H2O2 in the closed chemical system 
    with only dissocoation present remains constant

    """
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
           out_bin = '{"chem": {"rght": 1, "moms": ["O3_a", "H2O2_a", "SO2_a"], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0}}',)

    f = netcdf.netcdf_file(outfile,   "r")

    Mass      = [cm.M_SO2, cm.M_O3, cm.M_H2O2]
    chem_list = ["SO2",       "O3",    "H2O2"]
    eps       = [4e-9,       4e-11,      2e-4]
                                        #TODO why so big?

    # convert aqueous phase within droplets to gase phase (mole fraction)
    def aq_to_g(chem_aq, rhod, T, M, p):
        return chem_aq * rhod * cm.R * T / M  / p

    for chem in chem_list : 
        ini = f.variables[chem+"_g"][0] 
        #final gas phase = current gas phase + mass dissolved into droplets
        end = f.variables[chem+"_g"][-1] + \
              aq_to_g(f.variables[chem+"_a"][-1], f.variables["rhod"][-1], f.variables["T"][-1],\
                    Mass[chem_list.index(chem)], f.variables["p"][-1])

        assert np.isclose(end, ini, atol=0, rtol=eps[chem_list.index(chem)])
    
    # cleanup
    subprocess.call(["rm", outfile])
