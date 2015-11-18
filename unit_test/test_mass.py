# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")
from scipy.io import netcdf

import numpy as np
import pytest
import subprocess
import math

from parcel import parcel
from spectrum_plot import plot_spectrum
from libcloudphxx import common as cm
from functions import *

@pytest.fixture(scope="module")
def data(request):
    """
    Run parcel simulation and return opened netdcf file
    """
    outfile = "test_mass.nc"

    # Kreidenweis et al conditions
    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = rh_to_rv(RH_init, T_init, p_init)

    # STP conditions
    p_stp = 101325
    T_stp = 273.15 + 15

    rho_init = p_init / T_init / (r_init / (1.+r_init) * cm.R_v + 1./ (1.+r_init) * cm.R_d)
    rho_stp  = p_stp  / T_stp / cm.R_d
 
    # conversion from initial condidtion form Kreidenwieis et al paper 
    # to STP conditions needed for initialisation
    n_tot = 566.e6 * rho_stp / rho_init

    # running parcel model for open / closed chem system  ...
    parcel(dt = .1, sd_conc = 32168, outfreq = 1200, outfile = outfile, z_max = 600., w = .5,\
           T_0 = T_init, p_0 = p_init, r_0 = r_init, \
           chem_rho = 1.8e3, mean_r = .04e-6, gstdev = 2., n_tot = n_tot, out_bin = \
             '{"wradii": {"rght": 1e-4, "left": 1e-10, "drwt": "wet", "lnli": "lin", "nbin": 500, "moms": [0, 3]}, \
               "dradii": {"rght": 1e-4, "left": 1e-10, "drwt": "dry", "lnli": "lin", "nbin": 500, "moms": [0, 3]}}'
          )

    data = netcdf.netcdf_file(outfile, "r")

    # removing all netcdf files after all tests
    def removing_files():
        subprocess.call(["rm", outfile])

    #request.addfinalizer(removing_files)
    return data

def test_water_const(data, eps = 7e-15):
    """
    Check if the total water is preserved

    ini = water vapor mixing ratio at t = 0    + water for aerosol to reach equilibrum at t = 0
    end = water vapor mixing ratio at t = end  + water in all particles (cloud + aerosol) t = end

    """
    rho_w = cm.rho_w 
    rv    = data.variables["r_v"][:]
    mom3  = data.variables["wradii_m3"][:,:]
                                                    
    ini = mom3[0,:].sum()  * 4./3 * math.pi * rho_w + rv[0]
    end = mom3[-1,:].sum() * 4./3 * math.pi * rho_w + rv[-1]

#    assert np.isclose(end, ini, atol=0, rtol=eps), str((ini-end)/ini)

def test_dry_mass_const(data, eps = 1e-20):
    """
    Check if the total dry mass is preserved

    ini = dry particulate mass / kg dry air at t=0
    end = dry particulate mass / kg dry air at t=end

    """
    chem_rho = getattr(data, "chem_rho")
    mom3     = data.variables["dradii_m3"][:,:]
    rhod     = data.variables["rhod"][:]
    rv       = data.variables["r_v"][:]

    ini = mom3[0,:].sum()  * 4./3 * math.pi * chem_rho
    end = mom3[-1,:].sum() * 4./3 * math.pi * chem_rho

#    assert np.isclose(end, ini, atol=0, rtol=eps), str((ini-end)/ini)

    #epp = cm.R_d / cm.R_v
    rhod_parc_init = rhod[0]
    #rho_parc_init  = rhod[0] * (rv[0] * rv[0] + 1.) / (rv[0] + 1.)

    a = 2.375 * 1e-9 / rhod_parc_init

    print " "
    print " "

    print "n_tot = 566.e6, sd_conc=1024, r = 2.0444 "
    print "n_tot = 566.e6, sd_conc=2048, r = 2.0051 "
    print "n_tot = 566.e6, sd_conc=4096, r = 1.9359 "
    print "n_tot = 566.e6, sd_conc=8192, r = 1.9130 "
    print "n_tot = 566.e6, sd_conc=16384,r = 1.9145 "
    print "n_tot = 566.e6, sd_conc=32768,r = 1.9273 "
    print " "

    print "with correction: "
    print "n_tot ~ 601.e6, sd_conc=8192,  r = 2.0293 "
    print "n_tot ~ 601.e6, sd_conc=32768, r = 2.0459 "
    print "n_tot ~ 601.e6, sd_conc=65532, r = 2.0459 "
    print " "

    print "initial mixing ratio: ", ini, " [kg/kg dry air]"
    print "final mixing ratio:   ", end, " [kg/kg dry air]"
    print " "

    print "article init mixing ratio (/ rhod parc init):  ", a, " [kg / kg dry air]"
