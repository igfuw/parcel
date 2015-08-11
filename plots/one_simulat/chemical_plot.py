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

def plot_chem(data, output_folder = '', output_title = ''):

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

    plt.savefig(output_folder + output_title + "stat.svg")

    #-----------------------------------------------------------------------
    # plot non-dissociating chem species
    plt.figure(2, figsize=(18,14))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('O3   gas mole fraction [ppb]')  
    plots[1].set_xlabel('H2O2 gas mole fraction [ppt]')
    plots[2].set_xlabel('O3_a     kg / kg dry air')  
    plots[3].set_xlabel('H2O2_a   kg / kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t = f.variables["t"][spn_idx:]

      plots[0].plot(f.variables["O3_g"][spn_idx:]   * 1e9 , t, style[i])
      plots[0].legend(loc='upper right')
      plots[1].plot(f.variables["H2O2_g"][spn_idx:] * 1e12 , t, style[i])
      plots[2].plot(f.variables["O3_a"][spn_idx:]          , t, style[i])
      plots[3].plot(f.variables["H2O2_a"][spn_idx:]        , t, style[i])

    plt.savefig(output_folder + output_title + "O3_H2O2.svg")

    #-----------------------------------------------------------------------
    # plot pH , H+, OH-
    plt.figure(3, figsize=(12,7))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('average pH')
    plots[1].set_xlabel('H   kg / kg dry air')
    plots[2].set_xlabel('OH  kg / kg dry air')

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

    plt.savefig(output_folder + output_title + "pH.svg")

    #-----------------------------------------------------------------------
    # plot NH3_g, NH3_a, NH4+ 
    plt.figure(4, figsize=(12,7))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('NH3  gas mole fraction [ppt]')
    plots[1].set_xlabel('NH3_a  kg / kg dry air')
    plots[2].set_xlabel('NH4_a  kg / kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]
      plots[0].plot(f.variables["NH3_g"][spn_idx:]  * 1e12            , t, style[i])
      plots[1].plot(f.variables["NH3_a"][spn_idx:]                    , t, style[i])
      plots[2].plot(np.squeeze(f.variables["plt_ch_NH4_a"][spn_idx:]) , t, style[i])

    plt.savefig(output_folder + output_title  + "NH3.svg")
    #-----------------------------------------------------------------------
    # plot HNO3_g, NHO3_a, NO3+ 
    plt.figure(5, figsize=(12,7))
    plots = []

    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('HNO3  gas mole fraction [ppt]')
    plots[1].set_xlabel('HNO3_a  kg / kg dry air')
    plots[2].set_xlabel('NO3_a  kg / kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]
      plots[0].plot(f.variables["HNO3_g"][spn_idx:]  * 1e12           , t, style[i])
      plots[1].plot(f.variables["HNO3_a"][spn_idx:]                   , t, style[i])
      plots[2].plot(np.squeeze(f.variables["plt_ch_NO3_a"][spn_idx:]) , t, style[i])

    plt.savefig(output_folder + output_title  + "HNO3.svg")
 
    #-----------------------------------------------------------------------
    # plot CO2_g, CO2_a, HCO3+, CO3--
    plt.figure(6, figsize=(18,14))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('CO2  gas mole fraction [ppm]')
    plots[1].set_xlabel('CO2_a   kg / kg dry air')
    plots[2].set_xlabel('HCO3_a  kg / kg dry air')  
    plots[3].set_xlabel('CO3_a   kg / kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]
      plots[0].plot(f.variables["CO2_g"][spn_idx:]  * 1e6              , t, style[i])
      plots[1].plot(f.variables["CO2_a"][spn_idx:]                     , t, style[i])
      plots[2].plot(np.squeeze(f.variables["plt_ch_HCO3_a"][spn_idx:]) , t, style[i])
      plots[3].plot(np.squeeze(f.variables["plt_ch_CO3_a"][spn_idx:])  , t, style[i])

    plt.savefig(output_folder + output_title + "CO2.svg")
 
    #-----------------------------------------------------------------------
    # plot SO2_g, SO2_a, HSO3+, SO3--, HSO4-, SO4--, S_VI
    plt.figure(7, figsize=(18,14))
    plots = []

    for i in range(9):
      plots.append(plt.subplot(3,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('SO2    gas mole fraction [ppt]')
    plots[1].set_xlabel('SO2_a  kg / kg dry air')
    plots[3].set_xlabel('HSO3_a kg / kg dry air')  
    plots[4].set_xlabel('SO3_a  kg / kg dry air')
    plots[6].set_xlabel('HSO4_a kg / kg dry air')  
    plots[7].set_xlabel('SO4_a  kg / kg dry air')
    plots[8].set_xlabel('S_VI   kg / kg dry air')

    for ax in plots:
      ax.set_ylabel('t [s]')
    
    for i, f in f_out_chem.iteritems():
      t   = f.variables["t"][spn_idx:]
      plots[0].plot(f.variables["SO2_g"][spn_idx:]  * 1e12             , t, style[i])
      plots[1].plot(f.variables["SO2_a"][spn_idx:]                     , t, style[i])
      plots[3].plot(np.squeeze(f.variables["plt_ch_HSO3_a"][spn_idx:]) , t, style[i])
      plots[4].plot(np.squeeze(f.variables["plt_ch_SO3_a"][spn_idx:])  , t, style[i])
      plots[6].plot(np.squeeze(f.variables["plt_ch_HSO4_a"][spn_idx:]) , t, style[i])
      plots[7].plot(np.squeeze(f.variables["plt_ch_SO4_a"][spn_idx:])  , t, style[i])
      plots[8].plot(np.squeeze(f.variables["plt_ch_S_VI"][spn_idx:])   , t, style[i])

    if "chem_SO2_a" in f.variables.keys() and "radii_m0" in f.variables.keys():
        # convert aqueous phase S_IV within droplets to gase phase (mole fraction)
        S_IV = (np.sum(f.variables["chem_SO2_a"][spn_idx:] * f.variables["radii_m0"][spn_idx:],  axis=1) / common.M_SO2_H2O + \
                np.sum(f.variables["chem_HSO3_a"][spn_idx:] * f.variables["radii_m0"][spn_idx:], axis=1) / common.M_HSO3 + \
                np.sum(f.variables["chem_SO3_a"][spn_idx:] * f.variables["radii_m0"][spn_idx:],  axis=1) / common.M_SO3) \
           * f.variables["rhod"][spn_idx:] * common.R * f.variables["T"][spn_idx:] / f.variables["p"][spn_idx:]\
           + f.variables["SO2_g"][spn_idx:]

        plots[5].set_xlabel('S_IV   conc ppb')  
        plots[5].plot(S_IV * 1e9, t, style[i])

    if "radii_m3" in f.variables.keys() and "radii_m0" in f.variables.keys():
        plots[2].set_xlabel('LWC    g/kg')  
        plots[2].plot(np.sum(f.variables["radii_m3"][spn_idx:], axis=1) * 4. / 3 * math.pi * 998.2 * 1000, t, style[i])

    plt.savefig(output_folder + output_title + "SO2.svg")
    #-----------------------------------------------------------------------
 
def main():
    # initial condition
    SO2_g_init  = 200e-12 
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6 
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12
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
            SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
            chem_sys = 'open',   outfile="test_plot_chem_open.nc",\
            sd_conc = sd_conc,\
            chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
            out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
           chem_sys = 'closed', outfile="test_plot_chem_closed.nc",\
           sd_conc = sd_conc,\
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
           out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, outfile="test_plot_chem_off.nc",\
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
