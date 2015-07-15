import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf
import numpy as np
import pytest
import math
import subprocess

from libcloudphxx import common as cm
from parcel import parcel

def test_chem_closed(eps = 1e-11):
    """
    Checking if the total mass of SO_2, O_3 and H2O2 in the closed chemical system 
    with only dissocoation present remains constatnt

    """
    SO2_g_init  =  200e-12
    O3_g_init   =  50e-9
    H2O2_g_init =  500e-12
    outfreq = 20
    z_max = 200.
    outfile = "test_chem_closed.nc"

    parcel(dt = .1, z_max = z_max, outfreq = outfreq, SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            chem_sys = 'closed',   outfile = outfile,\
            chem_dsl = True, chem_dsc = False, chem_rct = False,\
            out_bin = ["chem:0/1/1/lin/wet/O3_a,H2O2_a,SO2_a"],)

    f = netcdf.netcdf_file(outfile,   "r")

    M_dry = cm.R / cm.R_d
    Mass      = [cm.M_SO2, cm.M_O3, cm.M_H2O2]
    chem_list = ["SO2",       "O3",    "H2O2"]
    for chem in chem_list : 
        all_chem  = (f.variables[chem+"_g"][0] - f.variables[chem+"_g"][:]) * \
                      Mass[chem_list.index(chem)] / M_dry - f.variables[chem+"_a"][:]
 
        for idx in all_chem:
            assert abs(idx) <= eps 

        #TODO - why this doesn't work?
        #assert np.isclose(abs(all_chem), 0, atol=0, rtol=eps).all() 

    subprocess.call(["rm", outfile])
