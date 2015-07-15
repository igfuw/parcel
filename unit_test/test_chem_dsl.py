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

def test_chem_dsl(eps = 1e-13):
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
            out_bin = ["radii:0/1/1/lin/wet/3", "chem:0/1/1/lin/wet/O3_a,H2O2_a,SO2_a"],)

    f = netcdf.netcdf_file(outfile,   "r")

    # average drop volume
    vol = np.squeeze(f.variables["radii_m3"][:]) * 4/3. * math.pi

    # check if:  dissolved chem spec     = partial pressure of chem_spec                  * Henry cnst * (convert to kg)
    diff_SO2  = f.variables["SO2_a"][:]  - f.variables["SO2_g"][:]  * f.variables["p"][:] * cm.H_SO2  * cm.M_SO2  * vol
    diff_O3   = f.variables["O3_a"][:]   - f.variables["O3_g"][:]   * f.variables["p"][:] * cm.H_O3   * cm.M_O3   * vol
    diff_H2O2 = f.variables["H2O2_a"][:] - f.variables["H2O2_g"][:] * f.variables["p"][:] * cm.H_H2O2 * cm.M_H2O2 * vol
 
    assert np.all(diff_SO2 < eps)
    assert np.all(diff_O3 < eps)
    assert np.all(diff_H2O2 < eps)
    
    # cleanup
    subprocess.call(["rm", outfile])
