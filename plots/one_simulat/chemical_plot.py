import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "../../")

from scipy.io import netcdf
import numpy as np
import pytest
import math
import subprocess

from parcel import parcel
from libcloudphxx import common
from functions import *

def plot_chem(data, output_folder = '', output_title = ''):

    def dn_calc(m_end, m_beg, M):
        return (m_end - m_beg) / M

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    f_out_chem = {}
    if "open" in data.iterkeys():
        f_out_chem["open"] = data["open"]
    if "closed" in data.iterkeys():
        f_out_chem["closed"] = data["closed"]

    style = {
      "open"   : "b.-",
      "closed" : "g.-",
      "off"    : "r.-"
    }
    spn_idx = 0
    #spn_idx = int(math.ceil(float(f_out_chem['open'].chem_spn)/float(f_out_chem['open'].outfreq)))

    #-----------------------------------------------------------------------
    # plot p, T, RH and dry/wet moments
    plt.figure(1, figsize=(18,14))
    plots = []

    for i in range(9):
      plots.append(plt.subplot(3,3,i+1))
                             #(rows, columns, number)
    plots[0].set_xlabel('p [hPa]')
    plots[1].set_xlabel('T [K]')
    plots[2].set_xlabel('RH')

    plots[3].set_xlabel('m0 1/kg dry air ')
    plots[4].set_xlabel('m1 mm/kg dry air ')
    plots[5].set_xlabel('lwc g/kg dry air')

    plots[6].set_xlabel('m0_dry 1/kg dry air ')
    plots[7].set_xlabel('m1_dry um/kg dry air ')
    plots[8].set_xlabel('mass cont ug/kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in data.iteritems():
      t = f.variables["t"][spn_idx:]

      plots[0].plot(f.variables["p"][spn_idx:] / 100.   , t, style[i], label=i)
      plots[0].legend(loc='upper right')
      plots[1].plot(f.variables["T"][spn_idx:]          , t, style[i])
      plots[2].plot(
	f.variables["RH"][spn_idx:]                     , t, style[i], 
	[f.variables["RH"][spn_idx:].max()] * t.shape[0], t, style[i]
      )
      plots[3].plot(np.squeeze(f.variables["plt_rw_m0"][spn_idx:]), t, style[i])
      plots[4].plot(\
          np.squeeze(f.variables["plt_rw_m1"][spn_idx:]) / np.squeeze(f.variables["plt_rw_m0"][spn_idx:]) * 1e3, t, style[i])
      plots[5].plot(\
          np.squeeze(f.variables["plt_rw_m3"][spn_idx:]) * 4. / 3 * math.pi * 998.2 * 1e3, t, style[i])
 
      plots[6].plot(np.squeeze(f.variables["plt_rd_m0"][spn_idx:]), t, style[i])
      plots[7].plot(\
          np.squeeze(f.variables["plt_rd_m1"][spn_idx:]) / np.squeeze(f.variables["plt_rd_m0"][spn_idx:]) * 1e6, t, style[i])
      plots[8].plot(\
          np.squeeze(f.variables["plt_rd_m3"][spn_idx:]) * 4./ 3 * math.pi * getattr(f, 'chem_rho')  * 1e9, t, style[i])

    plt.savefig(output_folder + output_title + "stat.pdf")

    #-----------------------------------------------------------------------
    # plot non-dissociating chem species
    plt.figure(2, figsize=(18,14))
    plots = []

    for i in range(6):
      plots.append(plt.subplot(2,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('dn O3 g') 
    plots[1].set_xlabel('dn O3 a')  
    plots[2].set_xlabel('dn O3 sum')
  
    plots[3].set_xlabel('dn H2O2 g')
    plots[4].set_xlabel('dn H2O2 a')
    plots[5].set_xlabel('dn H2O2 sum')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t = f.variables["t"][spn_idx:]

      dn_o3g   = dn_calc(f.variables["O3_g"][spn_idx:],   f.variables["O3_g"][spn_idx],   cm.M_O3)
      dn_o3a   = dn_calc(f.variables["O3_a"][spn_idx:],   f.variables["O3_a"][spn_idx],   cm.M_O3)
      dn_h2o2g = dn_calc(f.variables["H2O2_g"][spn_idx:], f.variables["H2O2_g"][spn_idx], cm.M_H2O2)
      dn_h2o2a = dn_calc(f.variables["H2O2_a"][spn_idx:], f.variables["H2O2_a"][spn_idx], cm.M_H2O2)

      plots[0].plot(dn_o3g, t, style[i])
      plots[0].legend(loc='upper right')
      plots[1].plot(dn_o3a, t, style[i])
      plots[2].plot(dn_o3g + dn_o3a, t, style[i])
      plots[3].plot(dn_h2o2g, t, style[i])
      plots[4].plot(dn_h2o2a, t, style[i])
      plots[5].plot(dn_h2o2g + dn_h2o2a, t, style[i])

    plt.savefig(output_folder + output_title + "O3_H2O2.pdf")

    #-----------------------------------------------------------------------
    # plot pH , H+, OH-
    plt.figure(3, figsize=(12,7))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('average pH')
    plots[1].set_xlabel('m H kg/kg dry air')
    plots[2].set_xlabel('m OH kg/kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]
      n_H = np.squeeze(f.variables["plt_ch_H"][spn_idx:]) / common.M_H
      vol = np.squeeze(f.variables["plt_rw_m3"][spn_idx:]) * 4/3. * math.pi * 1e3  #litres
      pH  = -1 * np.log10(n_H / vol)

      plots[0].plot(pH,                                              t, style[i])
      plots[1].plot(np.squeeze(f.variables["plt_ch_H"][spn_idx:])  , t, style[i])
      plots[2].plot(np.squeeze(f.variables["plt_ch_OH"][spn_idx:]) , t, style[i])

    plt.savefig(output_folder + output_title + "pH.pdf")

    #-----------------------------------------------------------------------
    # plot NH3_g, NH3_a, NH4+ 
    plt.figure(4, figsize=(12,7))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('dn NH3 g')
    plots[1].set_xlabel('dn NH3 a')
    plots[2].set_xlabel('dn NH4 a')
    plots[3].set_xlabel('sum of dn')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]

      dn_nh3g = dn_calc(f.variables["NH3_g"][spn_idx:], f.variables["NH3_g"][spn_idx], cm.M_NH3)
      dn_nh3a = dn_calc(f.variables["NH3_a"][spn_idx:], f.variables["NH3_a"][spn_idx], cm.M_NH3_H2O)
      dn_nh4  = dn_calc(np.squeeze(f.variables["plt_ch_NH4_a"][spn_idx:]),\
                        np.squeeze(f.variables["plt_ch_NH4_a"][spn_idx]), cm.M_NH4)

      plots[0].plot(dn_nh3g, t, style[i])
      plots[1].plot(dn_nh3a, t, style[i])
      plots[2].plot(dn_nh4, t, style[i])
      plots[3].plot(dn_nh3g + dn_nh3a + dn_nh4, t, style[i])

    plt.savefig(output_folder + output_title  + "NH3.pdf")
    #-----------------------------------------------------------------------
    # plot HNO3_g, NHO3_a, NO3+ 
    plt.figure(5, figsize=(12,7))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('dn HNO3 g')
    plots[1].set_xlabel('dn HNO3 a')
    plots[2].set_xlabel('dn NO3 a')  
    plots[3].set_xlabel('sum of dn')  

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]

      dn_hno3g = dn_calc(f.variables["HNO3_g"][spn_idx:], f.variables["HNO3_g"][spn_idx], cm.M_HNO3)
      dn_hno3a = dn_calc(f.variables["HNO3_a"][spn_idx:], f.variables["HNO3_a"][spn_idx], cm.M_HNO3)
      dn_no3   = dn_calc(np.squeeze(f.variables["plt_ch_NO3_a"][spn_idx:]), \
                         np.squeeze(f.variables["plt_ch_NO3_a"][spn_idx]), cm.M_NO3) 

      plots[0].plot(dn_hno3g, t, style[i])
      plots[1].plot(dn_hno3a, t, style[i])
      plots[2].plot(dn_no3,   t, style[i])
      plots[3].plot(dn_hno3g + dn_hno3a + dn_no3,   t, style[i])

    plt.savefig(output_folder + output_title  + "HNO3.pdf")
 
    #-----------------------------------------------------------------------
    # plot CO2_g, CO2_a, HCO3+, CO3--
    plt.figure(6, figsize=(18,14))
    plots = []

    for i in range(6):
      plots.append(plt.subplot(2,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('dn CO2 g')
    plots[1].set_xlabel('dn CO2 a')
    plots[2].set_xlabel('dn HCO3 a')  
    plots[3].set_xlabel('dn CO3 a')  
    plots[4].set_xlabel('sum of dn')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]

      dn_co2g = dn_calc(f.variables["CO2_g"][spn_idx:], f.variables["CO2_g"][spn_idx], cm.M_CO2)
      dn_co2a = dn_calc(f.variables["CO2_a"][spn_idx:], f.variables["CO2_a"][spn_idx], cm.M_CO2_H2O)
      dn_hco3 = dn_calc(np.squeeze(f.variables["plt_ch_HCO3_a"][spn_idx:]),\
                        np.squeeze(f.variables["plt_ch_HCO3_a"][spn_idx]), cm.M_HCO3)
      dn_co3  = dn_calc(np.squeeze(f.variables["plt_ch_CO3_a"][spn_idx:]),\
                        np.squeeze(f.variables["plt_ch_CO3_a"][spn_idx]), cm.M_CO3)

      plots[0].plot(dn_co2g, t, style[i])
      plots[1].plot(dn_co2a, t, style[i])
      plots[2].plot(dn_hco3, t, style[i])
      plots[3].plot(dn_co3,  t, style[i])
      plots[4].plot(dn_co2g + dn_co2a + dn_hco3 + dn_co3,  t, style[i])

    plt.savefig(output_folder + output_title + "CO2.pdf")
 
    #-----------------------------------------------------------------------
    # plot SO2_g, SO2_a, HSO3+, SO3--, HSO4-, SO4--, S_VI
    plt.figure(7, figsize=(18,14))
    plots = []

    for i in range(9):
      plots.append(plt.subplot(3,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('dn SO2 g [#/kg dry air]')
    plots[1].set_xlabel('dn SO2 a  #/kg dry air')
    plots[2].set_xlabel('gas vol.conc SO2 [ppb]')
    plots[2].set_xticks([0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    plots[2].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3'])
    plots[2].set_xlim([0., 0.3])

    plots[3].set_xlabel('dn HSO3_a # / kg dry air')  
    plots[4].set_xlabel('dn SO3_a  # / kg dry air')
    plots[5].set_xlabel('sum of dn')

    plots[6].set_xlabel('dn HSO4_a # / kg dry air')  
    plots[7].set_xlabel('dn SO4_a  # / kg dry air')
    plots[8].set_xlabel('dn S_VI   # / kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t    = f.variables["t"][spn_idx:]
      p    = f.variables["p"][spn_idx:]
      T    = f.variables["T"][spn_idx:]
      rhod = f.variables["rhod"][spn_idx:]
  
      dn_so2g = dn_calc(f.variables["SO2_g"][spn_idx:], f.variables["SO2_g"][spn_idx], cm.M_SO2)
      dn_so2a = dn_calc(f.variables["SO2_a"][spn_idx:], f.variables["SO2_a"][spn_idx], cm.M_SO2_H2O)
      dn_hso3 = dn_calc(np.squeeze(f.variables["plt_ch_HSO3_a"][spn_idx:]), \
                        np.squeeze(f.variables["plt_ch_HSO3_a"][spn_idx]), cm.M_HSO3)
      dn_so3  = dn_calc(np.squeeze(f.variables["plt_ch_SO3_a"][spn_idx:]),\
                        np.squeeze(f.variables["plt_ch_SO3_a"][spn_idx]), cm.M_SO3)
      dn_hso4 = dn_calc(np.squeeze(f.variables["plt_ch_HSO4_a"][spn_idx:]), \
                        np.squeeze(f.variables["plt_ch_HSO4_a"][spn_idx]), cm.M_HSO4)
      dn_so4  = dn_calc(np.squeeze(f.variables["plt_ch_SO4_a"][spn_idx:]),\
                        np.squeeze(f.variables["plt_ch_SO4_a"][spn_idx]), cm.M_SO4)
      dn_h2so4= dn_calc(np.squeeze(f.variables["plt_ch_S_VI"][spn_idx:]),\
                        np.squeeze(f.variables["plt_ch_S_VI"][spn_idx]), cm.M_H2SO4)

      plots[0].plot(dn_so2g, t, style[i])
      plots[1].plot(dn_so2a, t, style[i])
      plots[2].plot(mix_ratio_to_mole_frac(f.variables["SO2_g"][:], p, common.M_SO2, T, rhod) * 1e9, t)

      plots[3].plot(dn_hso3, t, style[i])
      plots[4].plot(dn_so3,  t, style[i])
      plots[5].plot(dn_so2g + dn_so2a + dn_hso3 + dn_so3 + dn_hso4 + dn_so4,  t, style[i])

      plots[6].plot(dn_hso4, t, style[i])
      plots[7].plot(dn_so4,  t, style[i])
      plots[8].plot(dn_h2so4,t, style[i])

#    if "chem_SO2_a" in f.variables.keys() and "radii_m0" in f.variables.keys():
#        S_IV = np.sum(f.variables["chem_SO2_a"][spn_idx:])  / cm.M_SO2_H2O + \
#               np.sum(f.variables["chem_HSO3_a"][spn_idx:]) / common.M_HSO3 + \
#               np.sum(f.variables["chem_SO3_a"][spn_idx:])  / common.M_SO3 + \
#               f.variables["SO2_g"][spn_idx:] / cm.M_SO2
#
#        plots[5].set_xlabel('n S_IV')  
#        plots[5].plot(S_IV, t, style[i])
#
#    if "radii_m3" in f.variables.keys() and "radii_m0" in f.variables.keys():
#        plots[2].set_xlabel('LWC    g/kg')  
#        plots[2].plot(np.sum(f.variables["radii_m3"][spn_idx:], axis=1) * 4. / 3 * math.pi * 998.2 * 1000, t, style[i])

    plt.savefig(output_folder + output_title + "SO2.pdf")
    #-----------------------------------------------------------------------
 
def main():

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

    # output
    z_max       = 200.
    dt          = .05
    outfreq     = int(z_max / dt / 100)
    w           = 1.
    sd_conc     = 2048.

    # turn on chemistry
    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    # define output for moments and chemistry
    out_bin_chem = '{"plt_rw":   {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_rd":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_ch":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1,\
                                  "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                           "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                           "CO2_a",  "HCO3_a", "CO3_a",\
                                           "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]\
                    }}'

    out_bin      = '{"plt_rw": {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_rd": {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]}}'

    # running parcel model for open / closed / off chem system
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
            T_0 = T_init, p_0 = p_init, r_0 = r_init,\
            SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
            chem_sys = 'open',   outfile="test_plot_chem_open.nc",\
            sd_conc = sd_conc,\
            chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
            out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
           chem_sys = 'closed', outfile="test_plot_chem_closed.nc",\
           sd_conc = sd_conc,\
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
           out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, outfile="test_plot_chem_off.nc",\
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g_0=0, O3_g_0=0, H2O2_g_0=0, out_bin = out_bin, sd_conc = sd_conc)

    # TODO - why do I have to repeat this import here?
    from scipy.io import netcdf
    data = {'open'   : netcdf.netcdf_file("test_plot_chem_open.nc",   "r"),\
            'closed' : netcdf.netcdf_file("test_plot_chem_closed.nc", "r"),\
            'off'    : netcdf.netcdf_file("test_plot_chem_off.nc",    "r")}

    plot_chem(data, output_folder = "../outputs", output_title = "/test_plot_chem_")

    for name, netcdf in data.iteritems():
        netcdf.close()
        subprocess.call(["rm", "test_plot_chem_" + name + ".nc"])

if __name__ == '__main__':
    main()
