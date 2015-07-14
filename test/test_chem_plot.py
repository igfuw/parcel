import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
import matplotlib.pyplot as plt
from scipy.io import netcdf
from parcel import parcel
import numpy as np
import pytest
import math

def test_chem_plot():

    SO2_g_init  = 200e-12 
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    outfreq = 1000    
    z_max = 250.
    dt = .01
    w  = 1.

    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    spn_idx = int(math.ceil(float(chem_spn)/float(outfreq)))
    #spn_idx = 2000
    #spn_idx = 0

    out_bin_chem = ["radii:0/1/1/lin/wet/0,1,3", 
                    "radiidry:0/1/1/lin/dry/0,1,3",
                    "chem:0/1/1/lin/wet/O3_a,H2O2_a,SO2_a,H,OH,HSO3_a,SO3_a,HSO4_a,SO4_a,S_VI"
                   ]

    out_bin  = ["radii:0/1/1/lin/wet/0,1,3", 
                "radiidry:0/1/1/lin/dry/0,1,3"
               ]

    # running parcel model for open / closed chem system  ...
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
            SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            chem_sys = 'open',   outfile="test_chem_open.nc",\
            chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
            out_bin = out_bin_chem)
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
             SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
             chem_sys = 'closed', outfile="test_chem_closed.nc",\
             chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
             out_bin = out_bin_chem)
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, SO2_g_0=0, O3_g_0=0, H2O2_g_0=0,\
             outfile="test_chem_off.nc", out_bin = out_bin)

    # ... plotting the results ...
    f_out = {
      "open"   : netcdf.netcdf_file("test_chem_open.nc",   "r"),
      "closed" : netcdf.netcdf_file("test_chem_closed.nc", "r"),
      "off"    : netcdf.netcdf_file("test_chem_off.nc",    "r")
    }
    f_out_chem = {
      "open"   : netcdf.netcdf_file("test_chem_open.nc",   "r"),
      "closed" : netcdf.netcdf_file("test_chem_closed.nc", "r")
    }

    style = {
      "open"   : "b.-",
      "closed" : "g.-",
      "off"    : "r.-"
    }

    plt.figure(1, figsize=(20,30))
    plots = []

    for i in range(24):
      plots.append(plt.subplot(8,3,i+1))

    plots[0].set_xlabel('p [hPa]')
    plots[1].set_xlabel('T [K]')
    plots[2].set_xlabel('RH')

    plots[3].set_xlabel('m0 1/kg dry air ')
    plots[4].set_xlabel('m1 m/kg dry air ')
    plots[5].set_xlabel('m3 m^3/kg dry air')

    plots[6].set_xlabel('m0_dry 1/kg dry air ')
    plots[7].set_xlabel('m1_dry m/kg dry air ')
    plots[8].set_xlabel('m3_dry m^3/kg dry air')

    plots[9].set_xlabel('average pH')
    plots[10].set_xlabel('H   kg / kg dry air')
    plots[11].set_xlabel('OH  kg / kg dry air')
 
    plots[12].set_xlabel('SO2  gas mole fraction [ppt]')
    plots[13].set_xlabel('O3   gas mole fraction [ppt]')  
    plots[14].set_xlabel('H2O2 gas mole fraction [ppb]')

    plots[15].set_xlabel('SO2_a    kg / kg dry air')
    plots[16].set_xlabel('O3_a     kg / kg dry air')  
    plots[17].set_xlabel('H2O2_a   kg / kg dry air')

    plots[18].set_xlabel('HSO3_a kg / kg dry air')  
    plots[19].set_xlabel('SO3_a kg / kg dry air')
    plots[20].set_xlabel('')  
 
    plots[21].set_xlabel('HSO4_a kg / kg dry air')  
    plots[22].set_xlabel('SO4_a  kg / kg dry air')
    plots[23].set_xlabel('S_VI   kg / kg dry air')

    for ax in plots:
      ax.set_ylabel('z [m]')

    for i, f in f_out.iteritems():
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
 
    for i, f in f_out_chem.iteritems():
      n_H = np.squeeze(f.variables["chem_H"][spn_idx:]) / common.M_H
      vol = np.squeeze(f.variables["radii_m3"][spn_idx:]) * 4/3. * math.pi * 1e3  #litres
      pH  = -1 * np.log10(n_H / vol)

      plots[9].plot(pH,  z, style[i])
      plots[10].plot(np.squeeze(f.variables["chem_H"][spn_idx:])  , z, style[i])
      plots[11].plot(np.squeeze(f.variables["chem_OH"][spn_idx:]) , z, style[i])

      plots[12].plot(f.variables["SO2_g"][spn_idx:]  * 1e12 , z, style[i])
      plots[13].plot(f.variables["O3_g"][spn_idx:]   * 1e12 , z, style[i])
      plots[14].plot(f.variables["H2O2_g"][spn_idx:] * 1e9  , z, style[i])
 
      plots[15].plot(f.variables["SO2_a"][spn_idx:]  , z, style[i])
      plots[16].plot(f.variables["O3_a"][spn_idx:]   , z, style[i])
      plots[17].plot(f.variables["H2O2_a"][spn_idx:] , z, style[i])

      plots[18].plot(np.squeeze(f.variables["chem_HSO3_a"][spn_idx:]) , z, style[i])
      plots[19].plot(np.squeeze(f.variables["chem_SO3_a"][spn_idx:])  , z, style[i])
  
      plots[21].plot(np.squeeze(f.variables["chem_HSO4_a"][spn_idx:]) , z, style[i])
      plots[22].plot(np.squeeze(f.variables["chem_SO4_a"][spn_idx:])  , z, style[i])
      plots[23].plot(np.squeeze(f.variables["chem_S_VI"][spn_idx:])   , z, style[i])

    plt.savefig("plot_chem.svg")
