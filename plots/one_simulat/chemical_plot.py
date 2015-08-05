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

def plot_chem(data, output_folder = "/outputs"):

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    # TODO - copy keys and values from data dictionary
    f_out_chem = {
      "open"   : netcdf.netcdf_file("test_chem_open.nc",   "r"),
      "closed" : netcdf.netcdf_file("test_chem_closed.nc", "r")
    }
    # TODO - copy keys from data dictionary
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
    plots[4].set_xlabel('m1 m/kg dry air ')
    plots[5].set_xlabel('m3 m^3/kg dry air')

    plots[6].set_xlabel('m0_dry 1/kg dry air ')
    plots[7].set_xlabel('m1_dry m/kg dry air ')
    plots[8].set_xlabel('m3_dry m^3/kg dry air')

    for ax in plots:
      ax.set_ylabel('z [m]')

    for i, f in data.iteritems():
      z = f.variables["z"][spn_idx:]
      plots[0].plot(f.variables["p"][spn_idx:] / 100.   , z, style[i], label=i)
      plots[0].legend(loc='upper right')
      plots[1].plot(f.variables["T"][spn_idx:]          , z, style[i])
      plots[2].plot(
	f.variables["RH"][spn_idx:]                     , z, style[i], 
	[f.variables["RH"][spn_idx:].max()] * z.shape[0], z, style[i]
      )
      plots[3].plot(np.squeeze(f.variables["radii_m0"][spn_idx:]), z, style[i])
      plots[4].plot(np.squeeze(f.variables["radii_m1"][spn_idx:]), z, style[i])
      plots[5].plot(np.squeeze(f.variables["radii_m3"][spn_idx:]), z, style[i])
 
      plots[6].plot(np.squeeze(f.variables["radiidry_m0"][spn_idx:]), z, style[i])
      plots[7].plot(np.squeeze(f.variables["radiidry_m1"][spn_idx:]), z, style[i])
      plots[8].plot(np.squeeze(f.variables["radiidry_m3"][spn_idx:]), z, style[i])

    plt.savefig(output_folder + "/plot_chem_stat.svg")

    #-----------------------------------------------------------------------
    # plot non-dissociating chem species
    plt.figure(2, figsize=(18,14))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('O3   gas mole fraction [ppt]')  
    plots[1].set_xlabel('H2O2 gas mole fraction [ppb]')
    plots[2].set_xlabel('O3_a     kg / kg dry air')  
    plots[3].set_xlabel('H2O2_a   kg / kg dry air')

    for ax in plots:
      ax.set_ylabel('z [m]')

    for i, f in f_out_chem.iteritems():
      z = f.variables["z"][spn_idx:]
      plots[0].plot(f.variables["O3_g"][spn_idx:]   * 1e12 , z, style[i])
      plots[0].legend(loc='upper right')
      plots[1].plot(f.variables["H2O2_g"][spn_idx:] * 1e9  , z, style[i])
      plots[2].plot(f.variables["O3_a"][spn_idx:]   , z, style[i])
      plots[3].plot(f.variables["H2O2_a"][spn_idx:] , z, style[i])

    plt.savefig(output_folder + "/plot_chem_O3_H2O2.svg")

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
      ax.set_ylabel('z [m]')

    for i, f in f_out_chem.iteritems():
      n_H = np.squeeze(f.variables["chem_H"][spn_idx:]) / common.M_H
      vol = np.squeeze(f.variables["radii_m3"][spn_idx:]) * 4/3. * math.pi * 1e3  #litres
      pH  = -1 * np.log10(n_H / vol)

      plots[0].plot(pH,  z, style[i])
      plots[1].plot(np.squeeze(f.variables["chem_H"][spn_idx:])  , z, style[i])
      plots[2].plot(np.squeeze(f.variables["chem_OH"][spn_idx:]) , z, style[i])

    plt.savefig(output_folder + "/plot_chem_pH.svg")

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
      ax.set_ylabel('z [m]')

    for i, f in f_out_chem.iteritems():
      plots[0].plot(f.variables["NH3_g"][spn_idx:]  * 1e12 , z, style[i])
      plots[1].plot(f.variables["NH3_a"][spn_idx:]  , z, style[i])
      plots[2].plot(np.squeeze(f.variables["chem_NH4_a"][spn_idx:]) , z, style[i])

    plt.savefig(output_folder + "/plot_chem_NH3.svg")
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
      ax.set_ylabel('z [m]')

    for i, f in f_out_chem.iteritems():
      plots[0].plot(f.variables["HNO3_g"][spn_idx:]  * 1e12 , z, style[i])
      plots[1].plot(f.variables["HNO3_a"][spn_idx:]  , z, style[i])
      plots[2].plot(np.squeeze(f.variables["chem_NO3_a"][spn_idx:]) , z, style[i])

    plt.savefig(output_folder + "/plot_chem_HNO3.svg")
 
    #-----------------------------------------------------------------------
    # plot CO2_g, CO2_a, HCO3+, CO3--
    plt.figure(6, figsize=(18,14))
    plots = []

    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('CO2  gas mole fraction [pptm]')
    plots[1].set_xlabel('CO2_a   kg / kg dry air')
    plots[2].set_xlabel('HCO3_a  kg / kg dry air')  
    plots[3].set_xlabel('CO3_a   kg / kg dry air')  
 
    for ax in plots:
      ax.set_ylabel('z [m]')

    for i, f in f_out_chem.iteritems():
      plots[0].plot(f.variables["CO2_g"][spn_idx:]  * 1e6 , z, style[i])
      plots[1].plot(f.variables["CO2_a"][spn_idx:]  , z, style[i])
      plots[2].plot(np.squeeze(f.variables["chem_HCO3_a"][spn_idx:]) , z, style[i])
      plots[3].plot(np.squeeze(f.variables["chem_CO3_a"][spn_idx:]) , z, style[i])

    plt.savefig(output_folder + "/plot_chem_CO2.svg")
 
    #-----------------------------------------------------------------------
    # plot SO2_g, SO2_a, HSO3+, SO3--, HSO4-, SO4--, S_VI
    plt.figure(7, figsize=(18,14))
    plots = []

    for i in range(9):
      plots.append(plt.subplot(3,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('SO2  gas mole fraction [ppt]')
    plots[1].set_xlabel('SO2_a    kg / kg dry air')
    plots[2].set_xlabel('')  
    plots[3].set_xlabel('HSO3_a kg / kg dry air')  
    plots[4].set_xlabel('SO3_a  kg / kg dry air')
    plots[5].set_xlabel('')  
    plots[6].set_xlabel('HSO4_a kg / kg dry air')  
    plots[7].set_xlabel('SO4_a  kg / kg dry air')
    plots[8].set_xlabel('S_VI   kg / kg dry air')

    for ax in plots:
      ax.set_ylabel('z [m]')

    for i, f in f_out_chem.iteritems():
      plots[0].plot(f.variables["SO2_g"][spn_idx:]  * 1e12 , z, style[i])
      plots[1].plot(f.variables["SO2_a"][spn_idx:]  , z, style[i])
      plots[3].plot(np.squeeze(f.variables["chem_HSO3_a"][spn_idx:]) , z, style[i])
      plots[4].plot(np.squeeze(f.variables["chem_SO3_a"][spn_idx:])  , z, style[i])
      plots[6].plot(np.squeeze(f.variables["chem_HSO4_a"][spn_idx:]) , z, style[i])
      plots[7].plot(np.squeeze(f.variables["chem_SO4_a"][spn_idx:])  , z, style[i])
      plots[8].plot(np.squeeze(f.variables["chem_S_VI"][spn_idx:])   , z, style[i])

    plt.savefig(output_folder + "/plot_chem_SO2.svg")

def main():
    # initial condition
    SO2_g_init  = 200e-12 
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6 
    NH3_g_init  = 0.   #100e-12
    HNO3_g_init = 0.   #100e-12
    outfreq     = 100  #100    
    z_max       = 200.
    dt          = .05
    w           = 1.

    # turn on chemistry
    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    # define output for moments and chemistry
    out_bin_chem = '{"radii":   {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                    "radiidry": {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                    "chem":     {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1,\
                                  "moms": ["O3_a", "H2O2_a",\
                                           "H", "OH",\
                                           "SO2_a", "HSO3_a", "SO3_a",\
                                           "HSO4_a", "SO4_a", "S_VI",\
                                           "CO2_a", "HCO3_a", "CO3_a",\
                                           "NH3_a", "NH4_a",\
                                           "HNO3_a", "NO3_a"]\
                    }}'
    out_bin      = '{"radii":   {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
                    "radiidry": {"rght": 1, "left": 0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]}}'

    # running parcel model for open / closed chem system  ...
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
            SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
            chem_sys = 'open',   outfile="test_chem_open.nc",\
            chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
            out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
             SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
             CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
             chem_sys = 'closed', outfile="test_chem_closed.nc",\
             chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
             out_bin = out_bin_chem)

    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, outfile="test_chem_off.nc",\
           SO2_g_0=0, O3_g_0=0, H2O2_g_0=0, out_bin = out_bin)

    # TODO - why do I have to repeat this import here?
    from scipy.io import netcdf
    data = {'open'   : netcdf.netcdf_file("test_chem_open.nc", "r"),\
            'closed' : netcdf.netcdf_file("test_chem_closed.nc", "r"),\
            'off'    : netcdf.netcdf_file("test_chem_off.nc", "r")}

    plot_chem(data, output_folder = "../outputs")

    for name, netcdf in data.iteritems():
        netcdf.close()
        subprocess.call(["rm", "test_chem_" + name + ".nc"])

if __name__ == '__main__':
    main()
