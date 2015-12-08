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

    plots[0].set_xlabel('n O3g 1/kg dry air') 
    plots[1].set_xlabel('n O3a 1/kg dry air')  
    plots[2].set_xlabel('n O3 - init 1/kg dry air')
  
    plots[3].set_xlabel('n H2O2g 1/kg dry air')
    plots[4].set_xlabel('n H2O2a 1/kg dry air')
    plots[5].set_xlabel('n H2O2 - init 1/kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t = f.variables["t"][spn_idx:]

      n_o3g   = f.variables["O3_g"][spn_idx:] / cm.M_O3
      n_o3a   = f.variables["O3_a"][spn_idx:] / cm.M_O3
      n_h2o2g = f.variables["H2O2_g"][spn_idx:] / cm.M_H2O2
      n_h2o2a = f.variables["H2O2_a"][spn_idx:] / cm.M_H2O2

      plots[0].plot(n_o3g, t, style[i])
      plots[0].legend(loc='upper right')
      plots[1].plot(n_o3a, t, style[i])
      plots[2].plot(n_o3g + n_o3a - n_o3g[0] - n_o3a[0], t, style[i])
      plots[3].plot(n_h2o2g, t, style[i])
      plots[4].plot(n_h2o2a, t, style[i])
      plots[5].plot(n_h2o2g + n_h2o2a - n_h2o2g[0] - n_h2o2a[0], t, style[i])

    plt.savefig(output_folder + output_title + "O3_H2O2.pdf")

    #-----------------------------------------------------------------------
    # plot pH , H+, OH-
    plt.figure(3, figsize=(12,7))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('average pH')
    plots[1].set_xlabel('n H  1/kg dry air')
    plots[2].set_xlabel('n OH 1/kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]
      n_H = np.squeeze(f.variables["plt_ch_H"][spn_idx:]) / common.M_H
      vol = np.squeeze(f.variables["plt_rw_m3"][spn_idx:]) * 4/3. * math.pi * 1e3  #litres
      pH  = -1 * np.log10(n_H / vol)

      plots[0].plot(pH,                                                        t, style[i])
      plots[1].plot(np.squeeze(f.variables["plt_ch_H"][spn_idx:])  / cm.M_H  , t, style[i])
      plots[2].plot(np.squeeze(f.variables["plt_ch_OH"][spn_idx:]) / cm.M_OH , t, style[i])

    plt.savefig(output_folder + output_title + "pH.pdf")

    #-----------------------------------------------------------------------
    # plot NH3_g, NH3_a, NH4+ 
    plt.figure(4, figsize=(12,7))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n NH3 g 1/kg dry air')
    plots[1].set_xlabel('n NH3 a 1/kg dry air')
    plots[2].set_xlabel('n NH4 a 1/kg dry air')
    plots[3].set_xlabel('n - init 1/kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]

      n_nh3g = f.variables["NH3_g"][spn_idx:] / cm.M_NH3
      n_nh3a = f.variables["NH3_a"][spn_idx:] / cm.M_NH3_H2O
      n_nh4  = np.squeeze(f.variables["plt_ch_NH4_a"][spn_idx:]) / cm.M_NH4

      plots[0].plot(n_nh3g, t, style[i])
      plots[1].plot(n_nh3a, t, style[i])
      plots[2].plot(n_nh4,  t, style[i])
      
      plots[3].plot(n_nh3g + n_nh3a + n_nh4 - n_nh3g[0] - n_nh3a[0] - n_nh4[0], t, style[i])

    plt.savefig(output_folder + output_title  + "NH3.pdf")

    #-----------------------------------------------------------------------
    # plot HNO3_g, NHO3_a, NO3+ 
    plt.figure(5, figsize=(12,7))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n HNO3 g 1/kg dry air')
    plots[1].set_xlabel('n HNO3 a 1/kg dry air')
    plots[2].set_xlabel('n NO3  a 1/kg dry air')  
    plots[3].set_xlabel('n - init 1/kg dry air')  

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]

      n_hno3g = f.variables["HNO3_g"][spn_idx:] / cm.M_HNO3
      n_hno3a = f.variables["HNO3_a"][spn_idx:] / cm.M_HNO3
      n_no3   = np.squeeze(f.variables["plt_ch_NO3_a"][spn_idx:]) / cm.M_NO3

      plots[0].plot(n_hno3g, t, style[i])
      plots[1].plot(n_hno3a, t, style[i])
      plots[2].plot(n_no3,   t, style[i])
      plots[3].plot(n_hno3g + n_hno3a + n_no3 - n_hno3g[0] - n_hno3a[0] - n_no3[0], t, style[i])

    plt.savefig(output_folder + output_title  + "HNO3.pdf")
 
    #-----------------------------------------------------------------------
    # plot CO2_g, CO2_a, HCO3+, CO3--
    plt.figure(6, figsize=(18,14))
    plots = []

    for i in range(6):
      plots.append(plt.subplot(2,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n CO2 g 1/kg dry air')
    plots[1].set_xlabel('n CO2 a 1/kg dry air')
    plots[2].set_xlabel('n HCO3 a 1/kg dry air')  
    plots[3].set_xlabel('n CO3 a 1/kg dry air')  
    plots[4].set_xlabel('n - int 1/kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]

      n_co2g = f.variables["CO2_g"][spn_idx:] / cm.M_CO2
      n_co2a = f.variables["CO2_a"][spn_idx:] / cm.M_CO2_H2O
      n_hco3 = np.squeeze(f.variables["plt_ch_HCO3_a"][spn_idx:]) / cm.M_HCO3
      n_co3  = np.squeeze(f.variables["plt_ch_CO3_a"][spn_idx:]) / cm.M_CO3

      plots[0].plot(n_co2g, t, style[i])
      plots[1].plot(n_co2a, t, style[i])
      plots[2].plot(n_hco3, t, style[i])
      plots[3].plot(n_co3,  t, style[i])
      plots[4].plot(n_co2g + n_co2a + n_hco3 + n_co3 - n_co2g[0] - n_co2a[0] - n_hco3[0] - n_co3[0],  t, style[i])

    plt.savefig(output_folder + output_title + "CO2.pdf")
 
    #-----------------------------------------------------------------------
    # plot SO2_g, SO2_a, HSO3+, SO3--, HSO4-, SO4--, S_VI
    plt.figure(7, figsize=(18,14))
    plots = []

    for i in range(9):
      plots.append(plt.subplot(3,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n SO2 g [1/kg dry air]')
    plots[1].set_xlabel('n SO2 a  1/kg dry air')
    plots[2].set_xlabel('gas vol.conc SO2 [ppb]')
    plots[2].set_xticks([0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    plots[2].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3'])
    plots[2].set_xlim([0., 0.3])

    plots[3].set_xlabel('n HSO3_a 1 / kg dry air')  
    plots[4].set_xlabel('n SO3_a  1 / kg dry air')
    plots[5].set_xlabel('n - init 1/kg dry air')

    plots[6].set_xlabel('n HSO4_a 1 / kg dry air')  
    plots[7].set_xlabel('n SO4_a  1 / kg dry air')
    plots[8].set_xlabel('n S_VI   1 / kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t    = f.variables["t"][spn_idx:]
      p    = f.variables["p"][spn_idx:]
      T    = f.variables["T"][spn_idx:]
      rhod = f.variables["rhod"][spn_idx:]
  
      n_so2g = f.variables["SO2_g"][spn_idx:] / cm.M_SO2
      n_so2a = f.variables["SO2_a"][spn_idx:] / cm.M_SO2_H2O
      n_hso3 = np.squeeze(f.variables["plt_ch_HSO3_a"][spn_idx:]) / cm.M_HSO3
      n_so3  = np.squeeze(f.variables["plt_ch_SO3_a"][spn_idx:])  / cm.M_SO3
      n_hso4 = np.squeeze(f.variables["plt_ch_HSO4_a"][spn_idx:]) / cm.M_HSO4
      n_so4  = np.squeeze(f.variables["plt_ch_SO4_a"][spn_idx:])  / cm.M_SO4
      n_h2so4= np.squeeze(f.variables["plt_ch_S_VI"][spn_idx:])   / cm.M_H2SO4

      plots[0].plot(n_so2g, t, style[i])
      plots[1].plot(n_so2a, t, style[i])
      plots[2].plot(mix_ratio_to_mole_frac(f.variables["SO2_g"][:], p, common.M_SO2, T, rhod) * 1e9, t)

      plots[3].plot(n_hso3, t, style[i])
      plots[4].plot(n_so3,  t, style[i])
      plots[5].plot(n_so2g + n_so2a + n_hso3 + n_so3 + n_hso4 + n_so4\
                  - n_so2g[0] + n_so2a[0] + n_hso3[0] + n_so3[0] + n_h2so4[0],  t, style[i])

      plots[6].plot(n_hso4, t, style[i])
      plots[7].plot(n_so4,  t, style[i])
      plots[8].plot(n_h2so4,t, style[i])

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
