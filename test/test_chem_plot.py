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
    O3_g_init   = 0.#50e-9
    H2O2_g_init = 500e-12
    outfreq = 10    
    z_max = 200.
    dt = .1
    w  = 1.

    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    spn_idx = int(math.ceil(float(chem_spn)/float(outfreq)))
    #spn_idx = 0

    out_wet = ["radii:0/1/1/lin/0,1,3", "chem:0/1/1/lin/O3_a,H2O2_a,SO2_a,H,OH,HSO3_a,SO3_a,HSO4_a,SO4_a,S_VI",]

    # running parcel model for open / closed chem system  ...
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            chem_sys = 'open',   outfile="test_chem_open.nc",\
            chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
            out_wet = out_wet)
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
             chem_sys = 'closed', outfile="test_chem_closed.nc",\
             chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
             out_wet = out_wet)
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq, SO2_g_0=0, O3_g_0=0, H2O2_g_0=0, outfile="test_chem_off.nc",\
             out_wet = ["radii:0/1/1/lin/0,1,3"])

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

    plt.figure(1, figsize=(20,24))
    plots = []

    for i in range(21):
      plots.append(plt.subplot(7,3,i+1))

    plots[0].set_xlabel('p [hPa]')
    plots[1].set_xlabel('T [K]')
    plots[2].set_xlabel('RH')

    plots[3].set_xlabel('m0 1/kg dry air ')
    plots[4].set_xlabel('m1 m/kg dry air ')
    plots[5].set_xlabel('m3 m^3/kg dry air')

    plots[6].set_xlabel('average pH')
    plots[7].set_xlabel('H   kg / kg dry air')
    plots[8].set_xlabel('OH  kg / kg dry air')
 
    plots[9].set_xlabel( 'SO2  gas mole fraction [ppt]')
    plots[10].set_xlabel('O3   gas mole fraction [ppt]')  
    plots[11].set_xlabel('H2O2 gas mole fraction [ppb]')

    plots[12].set_xlabel('SO2_a    kg / kg dry air')
    plots[13].set_xlabel('O3_a     kg / kg dry air')  
    plots[14].set_xlabel('H2O2_a   kg / kg dry air')

    plots[15].set_xlabel('HSO3_a kg / kg dry air')  
    plots[16].set_xlabel('SO3_a kg / kg dry air')
    plots[17].set_xlabel('')  
 
    plots[18].set_xlabel('HSO4_a kg / kg dry air')  
    plots[19].set_xlabel('SO4_a  kg / kg dry air')
    plots[20].set_xlabel('S_VI   kg / kg dry air')

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
 
    for i, f in f_out_chem.iteritems():
      n_H = np.squeeze(f.variables["chem_H"][spn_idx:]) / common.M_H
      vol = np.squeeze(f.variables["radii_m3"][spn_idx:]) * 4/3. * math.pi * 1e3  #litres
      pH  = -1 * np.log10(n_H / vol)
      plots[6].plot(pH,  z, style[i])
      plots[7].plot(np.squeeze(f.variables["chem_H"][spn_idx:])   , z, style[i])
      plots[8].plot(np.squeeze(f.variables["chem_OH"][spn_idx:]) , z, style[i])

      plots[9].plot(f.variables["SO2_g"][spn_idx:]  * 1e12 , z, style[i])
      plots[10].plot(f.variables["O3_g"][spn_idx:]   * 1e12 , z, style[i])
      plots[11].plot(f.variables["H2O2_g"][spn_idx:] * 1e9 , z, style[i])
 
      plots[12].plot(f.variables["SO2_a"][spn_idx:]   , z, style[i])
      plots[13].plot(f.variables["O3_a"][spn_idx:]   , z, style[i])
      plots[14].plot(f.variables["H2O2_a"][spn_idx:] , z, style[i])

      plots[15].plot(np.squeeze(f.variables["chem_HSO3_a"][spn_idx:])   , z, style[i])
      plots[16].plot(np.squeeze(f.variables["chem_SO3_a"][spn_idx:]) , z, style[i])
   
      plots[18].plot(np.squeeze(f.variables["chem_S_VI"][spn_idx:])   , z, style[i])
      plots[19].plot(np.squeeze(f.variables["chem_HSO4_a"][spn_idx:]) , z, style[i])
      plots[20].plot(np.squeeze(f.variables["chem_SO4_a"][spn_idx:]) , z, style[i])

    plt.savefig("plot_chem.svg")

