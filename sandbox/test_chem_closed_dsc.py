import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf
import numpy as np
import math
import subprocess

from parcel import parcel
from libcloudphxx import common as cm

def test_chem_closed():
    """
    Checking if the total mass of SO_2, O_3 and H2O2 in the closed chemical system 
    with only dissocoation present remains constant

    """
    #TODO - think of a similar test for reactions
    SO2_g_init  =  200e-12
    O3_g_init   =  50e-9
    H2O2_g_init =  500e-12
    outfreq     = 20
    z_max       = 200.
    outfile     = "test_chem_closed_dsc.nc"
    dt          = .01

    # run parcel
    parcel(dt = dt, z_max = z_max, outfreq = outfreq,\
           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
           chem_sys = 'closed',   outfile = outfile,\
           chem_dsl = True, chem_dsc = True, chem_rct = False,\
           out_bin = '{"chem": {"rght": 1, "left": 0, "drwt": "wet", "nbin": 1, "lnli": "lin",\
                       "moms": ["O3_a", "H2O2_a", "SO2_a", "HSO3_a", "SO3_a"]}}')

    f = netcdf.netcdf_file(outfile,   "r")

    Mass      = [cm.M_O3, cm.M_H2O2]
    chem_list = [   "O3",    "H2O2"]
    eps       = [  4e-11,      2e-4]

    # convert aqueous phase within droplets to gase phase (mole fraction)
    def aq_to_g(chem_aq, rhod, T, M, p):
        return chem_aq * rhod * cm.R * T / M  / p

    # check for O3 and H2O2 (they don't dissociate)
    for chem in chem_list : 
        ini = f.variables[chem+"_g"][0] 
        #final gas phase = current gas phase + mass dissolved into droplets
        end = f.variables[chem+"_g"][-1] + \
              aq_to_g(f.variables[chem+"_a"][-1], f.variables["rhod"][-1], f.variables["T"][-1],\
                    Mass[chem_list.index(chem)], f.variables["p"][-1])

        assert np.isclose(end, ini, atol=0, rtol=eps[chem_list.index(chem)])

    # check for SO2_g -> SO2_a + HSO3- + SO3--
    ini = f.variables["SO2_g"][0] 
    final_SO2_a = aq_to_g(f.variables["SO2_a"][-1],\
                          f.variables["rhod"][-1], \
                          f.variables["T"][-1], \
                          cm.M_SO2,\
                          f.variables["p"][-1])

    final_HSO3_a = aq_to_g(np.squeeze(f.variables["chem_HSO3_a"][-1]),\
                           f.variables["rhod"][-1],\
                           f.variables["T"][-1],\
                           cm.M_HSO3,\
                           f.variables["p"][-1])

    final_SO3_a = aq_to_g(np.squeeze(f.variables["chem_SO3_a"][-1]),\
                          f.variables["rhod"][-1],\
                          f.variables["T"][-1],\
                          cm.M_SO3,\
                          f.variables["p"][-1])

    #final gas phase = current gas phase + mass dissolved into droplets
    end = f.variables["SO2_g"][-1] + final_SO2_a + final_HSO3_a + final_SO3_a 

    print "  "
    print "M SO2 = " , cm.M_SO2
    print "M HSO3 = ", cm.M_HSO3
    print "M SO3 = " , cm.M_SO3
    print " "
    print "initial SO2 mole frac =   ", ini
    print "final SO2 mole frac =     ", f.variables["SO2_g"][-1]
    print "mole frac due to SO2_a =  ", final_SO2_a
    print "mole frac due to HSO3_a = ", final_HSO3_a
    print "mole frac due to SO3_a =  ", final_SO3_a
    print "sum at the end =          ", end
    print "sum - ini / ini =         ", (ini - end) / ini

    assert np.isclose(end, ini, atol=0, rtol=4e-9)
   
    # cleanup
    subprocess.call(["rm", outfile])
