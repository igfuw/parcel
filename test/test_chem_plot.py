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
    outfreq = 10
    z_max = 200.

    chem_dsl = True
    chem_dsc = True
    chem_rct = False

    out_wet = ["radii:0/1/1/lin/0,1,3", "chem:0/1/1/lin/O3_a,H2O2_a,SO2_a,H",]

    # running parcel model for open / closed chem system  ...
    parcel(dt = 1., z_max = z_max, outfreq = outfreq, SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            chem_sys = 'open',   outfile="test_chem_open.nc",\
            chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct,\
            out_wet = out_wet)
    parcel(dt = 1., z_max = z_max, outfreq = outfreq, SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
             chem_sys = 'closed', outfile="test_chem_closed.nc",\
             chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct,\
             out_wet = out_wet)
    parcel(dt = 1., z_max = z_max, outfreq = outfreq, SO2_g_0=0, O3_g_0=0, H2O2_g_0=0, outfile="test_chem_off.nc",\
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

    plt.figure(1, figsize=(18,15))
    plots = []

    for i in range(12):
      plots.append(plt.subplot(4,3,i+1))

    plots[0].set_xlabel('p [hPa]')
    plots[1].set_xlabel('T [K]')
    plots[2].set_xlabel('RH')

    plots[3].set_xlabel('H   mol / kg dry air')
    plots[4].set_xlabel('average volume l')
    plots[5].set_xlabel('average pH')
 
    plots[6].set_xlabel('SO2  gas mole fraction [ppt]')
    plots[7].set_xlabel('O3   gas mole fraction [ppt]')  
    plots[8].set_xlabel('H2O2 gas mole fraction [ppb]')

    plots[9].set_xlabel('SO2   kg / kg dry air')
    plots[10].set_xlabel('O3   kg / kg dry air')  
    plots[11].set_xlabel('H2O2 kg / kg dry air')

    for ax in plots:
      ax.set_ylabel('z [m]')

    for i, f in f_out.iteritems():
      z = f.variables["z"][:]
      plots[0].plot(f.variables["p"][:] / 100.   , z, style[i], label=i)
      plots[0].legend(loc='upper right')
      plots[1].plot(f.variables["T"][:]          , z, style[i])
      plots[2].plot(
	f.variables["RH"][:]                     , z, style[i], 
	[f.variables["RH"][:].max()] * z.shape[0], z, style[i]
      )
 
    for i, f in f_out_chem.iteritems():
      plots[6].plot(f.variables["SO2_g"][:]  * 1e12 , z, style[i])
      plots[7].plot(f.variables["O3_g"][:]   * 1e12 , z, style[i])
      plots[8].plot(f.variables["H2O2_g"][:] * 1e9 , z, style[i])
 
      plots[9].plot(f.variables["SO2_a"][:]   , z, style[i])
      plots[10].plot(f.variables["O3_a"][:]   , z, style[i])
      plots[11].plot(f.variables["H2O2_a"][:] , z, style[i])

      n_H = np.squeeze(f.variables["chem_H"][:]) / common.M_H
      vol = np.squeeze(f.variables["radii_m3"][:]) * 4/3. * math.pi * 1e3  #litres
      pH  = -1 * np.log10(n_H / vol)

      print i
      print pH     

      plots[3].plot(n_H, z, style[i])
      plots[4].plot(vol, z, style[i])
      plots[5].plot(pH,  z, style[i])

    plt.savefig("plot_chem.svg")

