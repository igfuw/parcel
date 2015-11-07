import sys
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
from chemical_plot import plot_chem
from functions import *

@pytest.fixture(scope="module")
def data(request):

    # initial condition
    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = rh_to_rv(RH_init, T_init, p_init)

    # calculate rhod for initial gas mixing ratio
    rhod_init   = rhod_calc(T_init, p_init, r_init)
    # initial condition for trace geses
    SO2_g_init  = mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
    O3_g_init   = mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
    H2O2_g_init = mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
    CO2_g_init  = mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
    NH3_g_init  = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
    HNO3_g_init = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

    # aerosol size distr.
    mean_r = .08e-6 / 2
    gstdev = 2./1
    n_tot  = 566.e6

    # chem process toggling
    chem_dsl = True
    chem_dsc = True
    chem_rct = False

    # output
    z_max       = 200
    dt          = .1
    w           = 1.
    outfreq     = int(z_max / dt / 25)
    sd_conc     = 4.
    outfile     = "test_chem_closed_dsc.nc"

    # run parcel
    parcel(dt = dt, z_max = z_max, outfreq = outfreq,\
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g = SO2_g_init,  O3_g = O3_g_init,   H2O2_g = H2O2_g_init,\
           CO2_g = CO2_g_init, NH3_g = NH3_g_init, HNO3_g = HNO3_g_init,\
           chem_sys = 'closed',   outfile = outfile,\
           mean_r = mean_r, gstdev = gstdev, n_tot = n_tot, sd_conc=sd_conc, \
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct,\
           out_bin = '{\
                  "chem"  : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500,\
                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                      "CO2_a",  "HCO3_a", "CO3_a",\
                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]},\
                  "radii" : {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500, "moms": [3]},\
                  "plt_rw": {"rght": 1,    "left": 0,    "drwt": "wet", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_rd": {"rght": 1,    "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_ch": {"rght": 1,    "left": 0,    "drwt": "dry", "lnli": "lin", "nbin": 1,\
                             "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                      "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                      "CO2_a",  "HCO3_a", "CO3_a",\
                                      "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'
    )

    data = netcdf.netcdf_file(outfile,   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", outfile])

    request.addfinalizer(removing_files)
    return data

def test_is_electroneutral(data, eps = 2e-7):
    """
    Check if after dissociation the electrical charge of cloud droplets is 0
    """
    # read the data
    # TODO - do it once for all the tests here
    m_H    = data.variables["chem_H"][-1, :]
    m_OH   = data.variables["chem_OH"][-1, :]
    m_HSO3 = data.variables["chem_HSO3_a"][-1, :]
    m_SO3  = data.variables["chem_SO3_a"][-1, :]
    m_HCO3 = data.variables["chem_HCO3_a"][-1, :]
    m_CO3  = data.variables["chem_CO3_a"][-1, :]
    m_NO3  = data.variables["chem_NO3_a"][-1, :]
    m_NH4  = data.variables["chem_NH4_a"][-1, :]
    m_HSO4 = data.variables["chem_HSO4_a"][-1, :]
    m_SO4  = data.variables["chem_SO4_a"][-1, :]

    # for each droplet check if positive ions == negative ions 
    for idx, val in np.ndenumerate(m_H):
      if m_H[idx] != 0 :
          # positive ions
          p1 = m_H[idx] / cm.M_H
          p2 = m_NH4[idx] / cm.M_NH4
          n_pos = p1 +p2
          # negative ions
          n1 = m_OH[idx] / cm.M_OH
          n2 = m_HSO3[idx] / cm.M_HSO3 + 2 * m_SO3[idx] / cm.M_SO3
          n3 = m_HCO3[idx] / cm.M_HCO3 + 2 * m_CO3[idx] / cm.M_CO3
          n4 = m_NO3[idx] / cm.M_NO3 
          n5 = m_HSO4[idx] / cm.M_HSO4 + 2 * m_SO4[idx] / cm.M_SO4
          n_neg = n1 + n2 + n3 + n4 + n5

          assert np.isclose(n_neg, n_pos, atol=0, rtol=eps),\
              "electoneutral assert error " + str((n_pos-n_neg)/n_pos) +\
              "\n n_pos = " + str(n_pos) + " and n_neg = " + str(n_neg) + \
              "\n with positive ions: " + \
              "\n     H   " + str(p1) + \
              "\n     NH4 "+ str(p2) + \
              "\n with negative ions: " + \
              "\n     OH  " + str(n1) + \
              "\n     S4  " + str(n2) + \
              "\n     CO2 " + str(n3) + \
              "\n     NO3 " + str(n4) + \
              "\n     S6  " + str(n5)

def test_is_mass_S6_const_with_dsl_dsc(data, eps=1e-15):
    """
    Check if the mas of H2SO4 remains constant.
    (There are no chemical reactions, so the mass should stay the same)
    """
    # read the data
    # TODO - do it once for all the tests here
    m_HSO4 = data.variables["chem_HSO4_a"][-1, :]
    m_SO4  = data.variables["chem_SO4_a"][-1, :]
    m_S6   = data.variables["chem_S_VI"][-1, :]

    m_S6_ini = data.variables["chem_S_VI"][0, :]

    # check mass
    n_S6_ini       = m_S6_ini.sum()  / cm.M_H2SO4
    n_S6_end       = m_S6.sum() / cm.M_H2SO4
    n_SO4_HSO4_end = m_SO4.sum() / cm.M_SO4 + m_HSO4.sum() / cm.M_HSO4

    assert np.isclose(n_S6_ini, n_S6_end, atol=0, rtol=eps), "1: " + str((n_S6_end - n_S6_ini) / n_S6_ini)
    assert np.isclose(n_S6_ini, n_SO4_HSO4_end, atol=0, rtol=eps), "2: " + str((n_SO4_HSO4_end - n_S6_ini) / n_S6_ini)

@pytest.mark.parametrize("ion", ["H2O", "SO2", "HSO3", "CO2", "HCO3", "NH3", "HNO3", "S6"])
def test_check_dissoc_constants(data, ion):
     """
     Check if the mass of chemical compounds agrees with the dissociation constants
     """
     # read the data
     # TODO - do it once for all the tests here
     V      = data.variables["radii_m3"][-1, :] * 4./3 * math.pi
     m_H    = data.variables["chem_H"][-1, :]
     m_OH   = data.variables["chem_OH"][-1, :]
 
     m_SO2  = data.variables["chem_SO2_a"][-1, :]
     m_HSO3 = data.variables["chem_HSO3_a"][-1, :]
     m_SO3  = data.variables["chem_SO3_a"][-1, :]
 
     m_HSO4 = data.variables["chem_HSO4_a"][-1, :]
     m_SO4  = data.variables["chem_SO4_a"][-1, :]
     m_S6   = data.variables["chem_S_VI"][-1, :]
 
     m_CO2  = data.variables["chem_CO2_a"][-1, :]
     m_HCO3 = data.variables["chem_HCO3_a"][-1, :]
     m_CO3  = data.variables["chem_CO3_a"][-1, :]

     m_HNO3 = data.variables["chem_HNO3_a"][-1, :]
     m_NO3  = data.variables["chem_NO3_a"][-1, :]
 
     m_NH4 = data.variables["chem_NH4_a"][-1, :]
     m_NH3 = data.variables["chem_NH3_a"][-1, :]

     def check_water(m_OH, m_H, vol, eps):
     # check for dissociation of water
     # dissociation constant of water K_H20 = [H][OH]
         max_error = 0
         for idx, vol in np.ndenumerate(V):
             if vol > 0:

                 H2O_dissoc    = m_OH[idx] / cm.M_OH / vol * m_H[idx] / cm.M_H / vol
                 current_error = abs((cm.K_H2O - H2O_dissoc) / H2O_dissoc)

                 if max_error < current_error:
                     max_error = current_error

         assert max_error <= eps, "K_H20 error = " + str(max_error)


     def check_ions(m_1, M_1, m_2, M_2, m_12, M_12, vol, teor_const, eps):
     # check for all species that have dissoc. constants K = [A][B]/[AB]
         max_error = 0
         for idx, vol in np.ndenumerate(V):
             if vol > 0:

                 dissoc_const  = (m_1[idx] / M_1 * m_2[idx] / M_2) / (m_12[idx] / M_12) / vol
                 current_error = abs((teor_const - dissoc_const) / teor_const)

                 if max_error < current_error:
                     max_error = current_error

         assert max_error <= eps, ion + " error = " + str(max_error)


     def check_S6_ions(m_HSO4, M_HSO4, m_SO4, M_SO4, m_S6, M_H2SO4, m_H, M_H, vol, K_HSO4, eps1, eps2):
     # check for dissociation of H2SO4 (assumes no non-dissociated H2SO4) 
         max_error_1 = 0
         max_error_2 = 0

         for idx, vol in np.ndenumerate(V):
             if vol > 0:

                 left_HSO4  = m_HSO4[idx] / M_HSO4 / vol
                 right_HSO4 = (m_H[idx] / M_H  / vol * m_S6[idx] / M_H2SO4 / vol) / (m_H[idx] / M_H / vol + K_HSO4)
                 left_SO4   = m_SO4[idx] / M_SO4 / vol
                 right_SO4  = (K_HSO4 * m_S6[idx] / M_H2SO4 / vol)  / (m_H[idx] / M_H / vol + K_HSO4)

                 current_error_1 = abs(left_HSO4 - right_HSO4) / left_HSO4
                 current_error_2 = abs(left_SO4 - right_SO4)   / left_SO4

                 if max_error_1 < current_error_1: max_error_1 = current_error_1
                 if max_error_2 < current_error_2: max_error_2 = current_error_2

         assert max_error_1 <= eps1, "K_HSO4 error = " + str(max_error_1)
         assert max_error_1 <= eps1, "K_SO4 error = " + str(max_error_2)
 

     if   ion == "H2O":  check_water(m_OH, m_H, V, 2e-2)

     elif ion == "SO2":  check_ions(m_H,   cm.M_H,   m_HSO3, cm.M_HSO3, m_SO2,  cm.M_SO2_H2O, V, cm.K_SO2,  2e-2)
     elif ion == "HSO3": check_ions(m_H,   cm.M_H,   m_SO3,  cm.M_SO3,  m_HSO3, cm.M_HSO3,    V, cm.K_HSO3, 3e-2) 
     elif ion == "CO2":  check_ions(m_H,   cm.M_H,   m_HCO3, cm.M_HCO3, m_CO2,  cm.M_CO2_H2O, V, cm.K_CO2,  2e-2)
     elif ion == "HCO3": check_ions(m_H,   cm.M_H,   m_CO3,  cm.M_CO3,  m_HCO3, cm.M_HCO3,    V, cm.K_HCO3, 3e-2)
     elif ion == "NH3":  check_ions(m_NH4, cm.M_NH4, m_OH,   cm.M_OH,   m_NH3,  cm.M_NH3_H2O, V, cm.K_NH3,  3e-2)
     elif ion == "HNO3": check_ions(m_H,   cm.M_H,   m_NO3,  cm.M_NO3,  m_HNO3, cm.M_HNO3,    V, cm.K_HNO3, 5e-4)
                                                                                                         # TODO - why so big?
                                                                                                         # check with > sd_conc
     elif ion == "S6":   check_S6_ions(m_HSO4, cm.M_HSO4, m_SO4, cm.M_SO4, m_S6, cm.M_H2SO4, m_H, cm.M_H,\
                                         V, cm.K_HSO4, 3e-2, 1e-3)
     else: assert False

@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2", "NH3", "HNO3"])
def test_is_mass_const_dsl_dsc(data, chem, eps =\
                                           {"SO2": 5e-4, "O3":9e-11, "H2O2": 4e-4, "CO2": 8e-10, "NH3": 2e-4, "HNO3":1e-4}):
                                            # TODO why so different?
     """
     Checking if the total number of moles in closed chemical system 
     with only dissocoation present remains constant
 
     """
     # check for O3 and H2O2 (they don't dissociate)
     if  chem in ["O3", "H2O2"] : 
         molar_mass = getattr(cm, "M_"+chem)
         ini = (data.variables[chem+"_g"][0] + data.variables[chem+"_a"][-1]) / molar_mass

         #final gas phase = current gas phase + mass dissolved into droplets
         end = (data.variables[chem+"_g"][-1] + data.variables[chem+"_a"][-1]) / molar_mass
 
     # check for NH3 -> NH4+ + OH-
     if chem == "NH3":
         ini = data.variables[chem+"_g"][0] / cm.M_NH3 +\
               data.variables["chem_"+chem+"_a"][0, :].sum()  / cm.M_NH3_H2O +\
               data.variables["chem_NH4_a"][0, :].sum() / cm.M_NH4

         #final gas phase = current gas phase + mass dissolved into droplets
         end = data.variables[chem+"_g"][-1] / cm.M_NH3 +\
               data.variables["chem_"+chem+"_a"][-1, :].sum() / cm.M_NH3_H2O +\
               data.variables["chem_NH4_a"][-1, :].sum() / cm.M_NH4
 
    # check for HNO3 -> H+ + NO3-
     if chem == "HNO3":
         ini = data.variables[chem+"_g"][0] / cm.M_HNO3 +\
               data.variables["chem_"+chem+"_a"][0, :].sum()  / cm.M_HNO3 +\
               data.variables["chem_NO3_a"][0, :].sum() / cm.M_NO3

         #final gas phase = current gas phase + mass dissolved into droplets
         end = data.variables[chem+"_g"][-1] / cm.M_HNO3 +\
               data.variables["chem_"+chem+"_a"][-1, :].sum() / cm.M_HNO3 +\
               data.variables["chem_NO3_a"][-1, :].sum() / cm.M_NO3
 
     # check for SO2_g -> SO2_a + HSO3- + SO3-- and the same for CO2
     if chem in ["SO2", "CO2"]:
         if chem == "SO2":
             M_gas  = cm.M_SO2
             M_aq   = cm.M_SO2_H2O
             M_ion1 = cm.M_HSO3
             M_ion2 = cm.M_SO3
         elif chem == "CO2":
             M_gas  = cm.M_CO2
             M_aq   = cm.M_CO2_H2O
             M_ion1 = cm.M_HCO3
             M_ion2 = cm.M_CO3
          
         ini = data.variables[chem+"_g"][0] / M_gas + \
               data.variables["chem_"+chem+"_a"][0, :].sum() / M_aq + \
               data.variables["chem_H"+chem.replace('2', '3')+"_a"][0, :].sum() / M_ion1 + \
               data.variables["chem_"+chem.replace('2','3')+"_a"][0, :].sum() / M_ion2

         #final gas phase = current gas phase + mass dissolved into droplets
         end = data.variables[chem+"_g"][-1] / M_gas + \
               data.variables["chem_"+chem+"_a"][-1, :].sum() / M_aq + \
               data.variables["chem_H"+chem.replace('2', '3')+"_a"][-1, :].sum() / M_ion1 + \
               data.variables["chem_"+chem.replace('2','3')+"_a"][-1, :].sum() / M_ion2

     assert np.isclose(end, ini, atol=0, rtol=eps[chem]), chem + " : " + str((ini-end)/ini)

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    data_to_plot = {'closed' : data}
    plot_chem(data_to_plot, output_folder="plots/outputs", output_title='/test_chem_closed_dsc_')

