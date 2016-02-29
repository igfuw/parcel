# This Python file uses the following encoding: utf-8
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
from libcloudphxx import common as cm
import functions as fn

def plot_chem(data, output_folder = '', output_title = ''):

    def dn_calc(m_end, m_beg, M):
        return (m_end - m_beg) / M

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    style = "b."     

    #spn_idx = 30
    spn_idx = 0

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
    plots[4].set_xlabel('m1 um/kg dry air ')
    plots[5].set_xlabel('lwc g/kg dry air')

    plots[6].set_xlabel('m0_dry 1/kg dry air ')
    plots[7].set_xlabel('m1_dry um/kg dry air ')
    plots[8].set_xlabel('mass cont ug/kg dry air')

    #plots[7].set_xlim([0.02, 0.04])
    #plots[7].set_xticks([0.02, 0.025, 0.03, 0.035, 0.04])

    for ax in plots:
      ax.set_ylabel('t [s]')

    t = data.variables["t"][spn_idx:]

    plots[0].plot(data.variables["p"][spn_idx:] / 100.   , t, style, label=i)
    plots[0].legend(loc='upper right')
    plots[1].plot(data.variables["T"][spn_idx:]          , t, style)
    plots[2].plot(
      data.variables["RH"][spn_idx:]                     , t, style, 
      [data.variables["RH"][spn_idx:].max()] * t.shape[0], t, style
    )
    plots[3].plot(np.squeeze(data.variables["plt_rw_m0"][spn_idx:]), t, style)
    plots[4].plot(\
        np.squeeze(data.variables["plt_rw_m1"][spn_idx:]) / np.squeeze(data.variables["plt_rw_m0"][spn_idx:]) * 1e6, t, style)
    plots[5].plot(\
        np.squeeze(data.variables["plt_rw_m3"][spn_idx:]) * 4. / 3 * math.pi * 998.2 * 1e3, t, style)
 
    plots[6].plot(np.squeeze(data.variables["plt_rd_m0"][spn_idx:]), t, style)
    plots[7].plot(\
        np.squeeze(data.variables["plt_rd_m1"][spn_idx:]) / np.squeeze(data.variables["plt_rd_m0"][spn_idx:]) * 1e6, t, style)
    plots[8].plot(\
        np.squeeze(data.variables["plt_rd_m3"][spn_idx:]) * 4./ 3 * math.pi * getattr(data, 'chem_rho')  * 1e9, t, style)

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

    t = data.variables["t"][spn_idx:]

    n_o3g   = data.variables["O3_g"][spn_idx:] / cm.M_O3
    n_o3a   = data.variables["O3_a"][spn_idx:] / cm.M_O3
    n_h2o2g = data.variables["H2O2_g"][spn_idx:] / cm.M_H2O2
    n_h2o2a = data.variables["H2O2_a"][spn_idx:] / cm.M_H2O2

    plots[0].plot(n_o3g, t, style)
    plots[0].legend(loc='upper right')
    plots[1].plot(n_o3a, t, style)
    plots[2].plot(n_o3g + n_o3a - n_o3g[0] - n_o3a[0], t, style)
    plots[3].plot(n_h2o2g, t, style)
    plots[4].plot(n_h2o2a, t, style)
    plots[5].plot(n_h2o2g + n_h2o2a - n_h2o2g[0] - n_h2o2a[0], t, style)

    plt.savefig(output_folder + output_title + "O3_H2O2.pdf")

    #-----------------------------------------------------------------------
    # plot pH , H+, OH-
    plt.figure(3, figsize=(12,8))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('average pH')
    plots[1].set_xlabel('n H  1/kg dry air')
    plots[2].set_xlabel('n OH 1/kg dry air')

    plots[2].ticklabel_format(axis = 'x', style='sci', scilimits=(-2,2))

    for ax in plots:
      ax.set_ylabel('t [s]')

    t   = data.variables["t"][spn_idx:]
    n_H = np.squeeze(data.variables["plt_ch_H"][spn_idx:]) / cm.M_H
    vol = np.squeeze(data.variables["plt_rw_m3"][spn_idx:]) * 4/3. * math.pi
    pH  = -1 * np.log10(n_H / vol / 1e3)
                                    # litres
    plots[0].plot(pH,                                                        t, style)
    plots[1].plot(np.squeeze(data.variables["plt_ch_H"][spn_idx:])  / cm.M_H  , t, style)
    plots[2].plot(cm.K_H2O / np.squeeze(data.variables["plt_ch_H"][spn_idx:]) * vol , t, style)

    plt.tight_layout()
    plt.savefig(output_folder + output_title + "pH.pdf")

    #-----------------------------------------------------------------------
    # plot NH3_g, NH3_a, NH4+ 
    plt.figure(4, figsize=(12,8))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n NH3 g 1/kg dry air')
    plots[1].set_xlabel('n N3 aq 1/kg dry air')
    plots[2].set_xlabel('n - init 1/kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    t    = data.variables["t"][spn_idx:]
    n_nh3g = data.variables["NH3_g"][spn_idx:] / cm.M_NH3
    n_N3 = np.squeeze(data.variables["plt_ch_NH3_a"][spn_idx:]) / cm.M_NH3_H2O

    plots[0].plot(n_nh3g, t, style)
    plots[1].plot(n_N3, t, style)
      
    plots[2].plot(n_nh3g + n_N3 - n_nh3g[0] - n_N3[0], t, style)

    plt.tight_layout()
    plt.savefig(output_folder + output_title  + "NH3.pdf")

    #-----------------------------------------------------------------------
    # plot HNO3_g, NHO3_a, NO3+ 
    plt.figure(5, figsize=(12,8))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n HNO3 g 1/kg dry air')
    plots[1].set_xlabel('n N5 aq  1/kg dry air')
    plots[2].set_xlabel('n - init 1/kg dry air')  

    for ax in plots:
      ax.set_ylabel('t [s]')

    t   = data.variables["t"][spn_idx:]
    n_hno3g = data.variables["HNO3_g"][spn_idx:] / cm.M_HNO3
    n_N5 = np.squeeze(data.variables["plt_ch_HNO3_a"][spn_idx:]) / cm.M_HNO3
 
    plots[0].plot(n_hno3g, t, style)
    plots[1].plot(n_N5,    t, style)
    plots[2].plot(n_hno3g + n_N5 - n_hno3g[0] - n_N5[0], t, style)

    plt.savefig(output_folder + output_title  + "HNO3.pdf")
 
    #-----------------------------------------------------------------------
    # plot CO2_g, CO2_a, HCO3+, CO3--
    plt.figure(6, figsize=(12,8))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n CO2 g 1/kg dry air')
    plots[1].set_xlabel('n C4 aq 1/kg dry air')
    plots[2].set_xlabel('n - int 1/kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    t   = data.variables["t"][spn_idx:]
    n_co2g = data.variables["CO2_g"][spn_idx:] / cm.M_CO2
    n_C4 = np.squeeze(data.variables["plt_ch_CO2_a"][spn_idx:]) / cm.M_CO2_H2O
 
    plots[0].plot(n_co2g, t, style)
    plots[1].plot(n_C4, t, style)
    plots[2].plot(n_co2g + n_C4 - n_co2g[0] - n_C4[0],  t, style)

    plt.savefig(output_folder + output_title + "CO2.pdf")
 
    #-----------------------------------------------------------------------
    # plot SO2_g, SO2_a, HSO3+, SO3--, HSO4-, SO4--, S_VI
    plt.figure(7, figsize=(18,14))
    plots = []

    for i in range(5):
      plots.append(plt.subplot(2,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('n SO2 g [1/kg dry air]')
    plots[1].set_xlabel('n S4 aq  1/kg dry air')

    plots[2].set_xlabel('gas vol.conc SO2 [ppb]')
    plots[2].set_xticks([0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    plots[2].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3'])
    plots[2].set_xlim([0., 0.3])

    plots[3].set_xlabel('n S6 [1 / kg dry air]')
    plots[4].set_xlabel('n - init 1/kg dry air')


    for ax in plots:
      ax.set_ylabel('t [s]')

    t    = data.variables["t"][spn_idx:]
    p    = data.variables["p"][spn_idx:]
    T    = data.variables["T"][spn_idx:]
    rhod = data.variables["rhod"][spn_idx:]
  
    n_so2g = data.variables["SO2_g"][spn_idx:] / cm.M_SO2
    n_S4 = np.squeeze(data.variables["plt_ch_SO2_a"][spn_idx:]) / cm.M_SO2_H2O
    n_S6 = np.squeeze(data.variables["plt_ch_S_VI"][spn_idx:])  / cm.M_H2SO4

    plots[0].plot(n_so2g, t, style)
    plots[1].plot(n_S4, t, style)
    plots[2].plot(fn.mix_ratio_to_mole_frac(data.variables["SO2_g"][:], p, cm.M_SO2, T, rhod) * 1e9, t)
    plots[3].plot(n_S6, t, style)
    plots[4].plot(n_so2g + n_S4 + n_S6 - n_so2g[0] - n_S4[0] - n_S6[0],  t, style)

    plt.savefig(output_folder + output_title + "SO2.pdf")
    #-----------------------------------------------------------------------
 
def main():

    # initial condition
    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = fn.rh_to_rv(RH_init, T_init, p_init) 

    # calculate rhod for initial gas mixing ratio
    rhod_init   = fn.rhod_calc(T_init, p_init, r_init)
    # initial condition for trace geses
    SO2_g_init  = fn.mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
    O3_g_init   = fn.mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
    H2O2_g_init = fn.mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
    CO2_g_init  = fn.mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
    NH3_g_init  = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
    HNO3_g_init = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

    # output
    z_max       = 20. #200.
    dt          = .05
    outfreq     = int(z_max / dt / 100)
    w           = 1.
    sd_conc     = 2048

    # turn on chemistry
    chem_dsl = True
    chem_dsc = True
    chem_rct = True

    # define output for moments and chemistry
    out_bin_chem = '{"plt_rw":   {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_rd":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_ch":   {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1,\
                                  "moms": ["O3_a",   "H2O2_a", "H", "SO2_a", "S_VI", "CO2_a",  "NH3_a",  "HNO3_a"]}}'

    out_bin      = '{"plt_rw": {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                     "plt_rd": {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]}}'

    # run parcel model
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g = SO2_g_init, O3_g  = O3_g_init,  H2O2_g = H2O2_g_init,\
           CO2_g = CO2_g_init, NH3_g = NH3_g_init, HNO3_g = HNO3_g_init,\
           outfile="test_plot_chem_closed.nc",\
           sd_conc = sd_conc,\
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, \
           out_bin = out_bin_chem)

    data = netcdf.netcdf_file("test_plot_chem.nc", "r")

    plot_chem(data, output_folder = "../outputs", output_title = "/test_plot_chem_")

    netcdf.close()
    subprocess.call(["rm", "test_plot_chem.nc"])

if __name__ == '__main__':
    main()
