import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf
import numpy as np
import math
import subprocess

from libcloudphxx import common as cm
from parcel import parcel
                        #TODO is it small enough?
def test_chem_dsl(eps = 2e-5):
    """
    Checking if dissolving chemical compounds into cloud droplets follows Henrys law
    http://www.henrys-law.org/

    """
    SO2_g_init  =  200e-12 
    O3_g_init   =  50e-9
    H2O2_g_init =  500e-12
    outfreq = 20
    z_max = 200.
    outfile = "test_chem_dsl.nc"

    # running parcel model for open chem system  and only for dissolving chem species into droplets
    parcel(dt = .5, z_max = z_max, outfreq = outfreq, SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            chem_sys = 'open',   outfile = outfile,\
            chem_dsl = True, chem_dsc = False, chem_rct = False,\
            out_bin = '{"radii": {"rght": 1, "moms": [3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0}, "chem": {"rght": 1, "moms": ["O3_a", "H2O2_a", "SO2_a"], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0}}')

    f = netcdf.netcdf_file(outfile,   "r")

    # average drop volume
    vol = np.squeeze(f.variables["radii_m3"][:]) * 4/3. * math.pi

    def henry_checker(conc_aq, conc_g, p, H, M, vol):
    # dissolved  = partial prs * Henry_const * molar mass * drop volume
        henry_aq = conc_g[1:] * p[1:] * H * M * vol[1:]
        conc_aq  = conc_aq[1:]

        assert np.isclose(conc_aq, henry_aq, atol=0, rtol=eps).all()
 
    henry_checker(f.variables["SO2_a"][:],  f.variables["SO2_g"][:],  f.variables["p"][:], cm.H_SO2,  cm.M_SO2, vol)
    henry_checker(f.variables["O3_a"][:],   f.variables["O3_g"][:],   f.variables["p"][:], cm.H_O3,   cm.M_O3, vol)
    henry_checker(f.variables["H2O2_a"][:], f.variables["H2O2_g"][:], f.variables["p"][:], cm.H_H2O2, cm.M_H2O2, vol)
   
    # cleanup
    subprocess.call(["rm", outfile])
