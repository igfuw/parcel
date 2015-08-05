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

    #TODO - think of a similar test for reactions
    SO2_g_init  = 200e-12
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6
    NH3_g_init  = 0.      # 100e-12
    HNO3_g_init = 0.      # 100e-12 #TODO - add HNO3!!!

    outfreq     = 1000
    z_max       = 100.
    outfile     = "test_chem_closed_dsc.nc"
    dt          = .01

    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = cm.eps * RH_init * cm.p_vs(T_init) / (p_init - RH_init * cm.p_vs(T_init))

    mean_r = .08e-6 / 2
    gstdev = 2./1
    n_tot  = 566.e6

    # run parcel
    parcel(dt = dt, z_max = z_max, outfreq = outfreq,\
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g_0 = SO2_g_init,  O3_g_0 = O3_g_init,   H2O2_g_0 = H2O2_g_init,\
           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
           chem_sys = 'closed',   outfile = outfile,\
           mean_r = mean_r, gstdev = gstdev, n_tot = n_tot,\
           chem_dsl = True, chem_dsc = True, chem_rct = False,\
           out_bin = '{"chem": {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "nbin": 500, "lnli": "log",\
                                "moms": ["O3_a", "H2O2_a", "SO2_a",\
                                         "CO2_a", "HNO3_a", "NH3_a",\
                                         "H", "OH",\
                                         "HSO3_a", "SO3_a",\
                                         "HCO3_a", "CO3_a",\
                                         "NH4_a", "NO3_a",\
                                         "HSO4_a", "SO4_a", "S_VI"\
                                        ]},\
                       "radii":{"rght": 1e-4, "left": 1e-9, "drwt": "wet", "nbin": 500, "lnli": "log",\
                                "moms": [3]}}')

    data = netcdf.netcdf_file(outfile,   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", outfile])

    #request.addfinalizer(removing_files)
    return data

def test_is_electroneutral(data, eps = 9e-10):
    """
    Check if after dissociation the electrical charge of cloud droplets is 0
    """
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

    # positive ions
    n_pos = m_H.sum() / cm.M_H + m_NH4.sum() / cm.M_NH4
    # negative ions
    n_neg = m_OH.sum() / cm.M_OH +\
            m_HSO3.sum() / cm.M_HSO3 + 2 * m_SO3.sum() / cm.M_SO3 +\
            m_HCO3.sum() / cm.M_HCO3 + 2 * m_CO3.sum() / cm.M_CO3 +\
            m_NO3.sum() / cm.M_NO3 +\
            m_HSO4.sum() / cm.M_HSO4 + 2 * m_SO4.sum() / cm.M_SO4

    assert np.isclose(n_neg, n_pos, atol=0, rtol=eps), str((n_pos-n_neg)/n_pos)

def test_is_mass_S6_const_with_dsl_dsc(data, eps=1e-10):
    """
    Check if the mas of H2SO4 remains constant.
    (The initial aerosol is assumed to be made of H2SO4. There are no chem. reactions, so the mass should stay the same)
    """
    m_HSO4 = data.variables["chem_HSO4_a"][-1, :]
    m_SO4  = data.variables["chem_SO4_a"][-1, :]
    m_S6   = data.variables["chem_S_VI"][-1, :]

    m_S6_ini = data.variables["chem_S_VI"][0, :]

    n_S6_ini       = m_S6_ini.sum()  / cm.M_H2SO4
    n_S6_end       = m_S6.sum() / cm.M_H2SO4
    n_SO4_HSO4_end = m_SO4.sum() / cm.M_SO4 + m_HSO4.sum() / cm.M_HSO4

    assert np.isclose(n_S6_ini, n_S6_end, atol=0, rtol=eps)
    assert np.isclose(n_S6_ini, n_SO4_HSO4_end, atol=0, rtol=eps)

def test_check_dissoc_constants(data):
     """
     Check if the mass of chemical compounds agrees with the dissociation constants
     """
     #TODO - cleanup
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
 
     m_NH4  = data.variables["chem_NH4_a"][-1, :]
     m_NH3 = data.variables["chem_NH3_a"][-1, :]
 
     global num
 
     num = 0
 
     H2O_dissoc = 0
     SO2_H2O_dissoc = 0
     HSO3_dissoc = 0
     CO2_H2O_dissoc = 0
     HCO3_dissoc = 0
     NH3_H2O_dissoc = 0
 
     left_HSO4  = 0
     right_HSO4 = 0
     left_SO4 = 0
     right_SO4 = 0
 
     for idx, vol in np.ndenumerate(V):
         if vol > 0:
             num += 1
             H2O_dissoc += m_OH[idx] / cm.M_OH / vol * m_H[idx] / cm.M_H / vol
             SO2_H2O_dissoc += (m_H[idx] / cm.M_H * m_HSO3[idx] / cm.M_HSO3) / (m_SO2[idx] / cm.M_SO2_H2O) / vol
             HSO3_dissoc    += (m_H[idx] / cm.M_H * m_SO3[idx]  / cm.M_SO3)  / (m_HSO3[idx] / cm.M_HSO3) / vol
             CO2_H2O_dissoc += (m_H[idx] / cm.M_H * m_HCO3[idx] / cm.M_HCO3) / (m_CO2[idx] / cm.M_CO2_H2O) / vol
             HCO3_dissoc    += (m_H[idx] / cm.M_H * m_CO3[idx]  / cm.M_CO3)  / (m_HCO3[idx] / cm.M_HCO3) / vol
#             NH3_H2O_dissoc += (m_NH4[idx] / cm.M_NH4 * m_OH[idx] / cm.M_OH) / (m_NH3[idx] / cm.M_NH3_H2O) / vol
 
             left_HSO4  += m_HSO4[idx] / cm.M_HSO4 / vol
             right_HSO4 += (m_H[idx] / cm.M_H  / vol * m_S6[idx] / cm.M_H2SO4 / vol) / (m_H[idx] / cm.M_H / vol + cm.K_HSO4)
             left_SO4  += m_SO4[idx] / cm.M_SO4 / vol
             right_SO4 += (cm.K_HSO4 * m_S6[idx] / cm.M_H2SO4 / vol)  / (m_H[idx] / cm.M_H / vol + cm.K_HSO4)
 
     assert np.isclose(cm.K_H2O, H2O_dissoc/num, atol=0, rtol=5e-4),       str((cm.K_H2O - H2O_dissoc/num) / (H2O_dissoc/num))
     assert np.isclose(cm.K_SO2, SO2_H2O_dissoc/num , atol=0, rtol=5e-4),  str((cm.K_SO2 - SO2_H2O_dissoc/num) /cm.K_SO2 )
     assert np.isclose(cm.K_HSO3, HSO3_dissoc/num, atol=0, rtol=9e-4),     str((cm.K_HSO3 - HSO3_dissoc/num) / cm.K_HSO3)
 
     assert np.isclose(cm.K_CO2, CO2_H2O_dissoc/num , atol=0, rtol=3e-4),     str((cm.K_CO2 - CO2_H2O_dissoc/num) /cm.K_CO2 )
     assert np.isclose(cm.K_HCO3, HCO3_dissoc/num, atol=0, rtol=6e-4),        str((cm.K_HCO3 - HCO3_dissoc/num) / cm.K_HCO3)
#     assert np.isclose(cm.K_NH3, NH3_H2O_dissoc/num, atol=0, rtol=3e-4),      str((cm.K_NH3 - NH3_H2O_dissoc/num) / cm.K_NH3)
 
     assert np.isclose(left_HSO4/num, right_HSO4/num, atol=0, rtol=4e-10), str((left_HSO4/num - right_HSO4/num) / (left_HSO4/num))
     assert np.isclose(left_SO4/num, right_SO4/num, atol=0, rtol=3e-8),    str((left_SO4/num - right_SO4/num) / (left_SO4/num))
 
@pytest.mark.parametrize("chem", ["SO2", "O3", "H2O2", "CO2"]) #, "NH3"])
def test_is_mass_const_dsl_dsc(data, chem, eps = {"SO2": 4e-5, "O3":4e-11, "H2O2": 2e-4, "CO2": 2e-10, "NH3": 2e-4}):   
                                                   # TODO why so different?
     """
     Checking if the total mass of SO_2, O_3 and H2O2 in the closed chemical system 
     with only dissocoation present remains constant
 
     """
     if chem in ["SO2", "CO2", "NH3"]:
         molar_mass = getattr(cm, "M_"+chem+"_H2O")
     else: molar_mass = getattr(cm, "M_"+chem)
 
    # convert aqueous phase within droplets to gase phase (mole fraction)
     def aq_to_g(chem_aq, rhod, T, M, p):
         return chem_aq * rhod * cm.R * T / M  / p
 
     # check for O3 and H2O2 (they don't dissociate)
     if  chem in ["O3", "H2O2"] : 
         ini = data.variables[chem+"_g"][0] 
         #final gas phase = current gas phase + mass dissolved into droplets
         end = data.variables[chem+"_g"][-1] + \
               aq_to_g(data.variables[chem+"_a"][-1], data.variables["rhod"][-1], data.variables["T"][-1],\
                     molar_mass, data.variables["p"][-1])
 
     # check for NH3 -> NH4+ + OH-
     if chem == "NH3":
         ini = data.variables[chem+"_g"][0]
         #final gas phase = current gas phase + mass dissolved into droplets
         end = data.variables[chem+"_g"][-1] +\
               aq_to_g(
                   data.variables["chem_"+chem+"_a"][-1, :].sum(),
                   data.variables["rhod"][-1],
                   data.variables["T"][-1],
                   cm.M_NH3_H2O,
                   data.variables["p"][-1]
               )+\
               aq_to_g(
                   data.variables["chem_NH4_a"][-1, :].sum(),
                   data.variables["rhod"][-1],
                   data.variables["T"][-1],
                   cm.M_NH4,
                   data.variables["p"][-1]
               )
 
     # check for SO2_g -> SO2_a + HSO3- + SO3-- and the same for CO2
     if chem in ["SO2", "CO2"]:
         if chem == "SO2":
             M_ion1 = cm.M_HSO3
             M_ion2 = cm.M_SO3
         elif chem == "CO2":
             M_ion1 = cm.M_HCO3
             M_ion2 = cm.M_CO3
          
         ini = data.variables[chem+"_g"][0]
         #final gas phase = current gas phase + mass dissolved into droplets
         end = data.variables[chem+"_g"][-1] + \
               aq_to_g(
                   data.variables["chem_"+chem+"_a"][-1, :].sum(),
                   data.variables["rhod"][-1], 
                   data.variables["T"][-1], 
                   molar_mass,
                   data.variables["p"][-1]
               ) + \
               aq_to_g(
                   data.variables["chem_H"+chem.replace('2', '3')+"_a"][-1, :].sum(),
                   data.variables["rhod"][-1],
                   data.variables["T"][-1],
                   M_ion1,
                   data.variables["p"][-1]
               ) + \
               aq_to_g(
                   data.variables["chem_"+chem.replace('2','3')+"_a"][-1, :].sum(),
                   data.variables["rhod"][-1],
                   data.variables["T"][-1],
                   M_ion2,
                   data.variables["p"][-1])

 
     assert np.isclose(end, ini, atol=0, rtol=eps[chem]), chem + " : " + str((ini-end)/ini)
