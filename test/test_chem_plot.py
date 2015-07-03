import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
import matplotlib.pyplot as plt
from scipy.io import netcdf
from parcel import parcel
import numpy as np
import pytest

def test_chem_plot():
    # running parcel model for open / closed chem system  ...
#    parcel(dt = 1., outfreq = 10, SO2_g_0 = 200e-12, O3_g_0 = 50e-9, H2O2_g_0 = 500e-12,\
#            chem_sys = 'open',   outfile="test_chem_open.nc")
    parcel(dt = 1., outfreq = 10, SO2_g_0 = 200e-12, O3_g_0 = 50e-9, H2O2_g_0 = 500e-12,\
             chem_sys = 'closed', outfile="test_chem_closed.nc")
    parcel(dt = 1., outfreq = 10, SO2_g_0 = 200e-12, O3_g_0 = 50e-9, H2O2_g_0 = 500e-12,\
            chem_sys = 'none', outfile="test_chem_none.nc")
    parcel(dt = 1., outfreq = 10, SO2_g_0=0, O3_g_0=0, H2O2_g_0=0, outfile="test_chem_off.nc")

    # ... plotting the results ...
    f_out = {
#      "open"   : netcdf.netcdf_file("test_chem_open.nc",   "r"),
      "closed" : netcdf.netcdf_file("test_chem_closed.nc", "r"),
      "none"   : netcdf.netcdf_file("test_chem_none.nc", "r"),
      "off"    : netcdf.netcdf_file("test_chem_off.nc", "r")
    }

    style = {
#      "open"   : "b.-",
      "closed" : "g.-",
      "none"   : "r.-",
      "off"    : "m.-"
    }

    plt.figure(1, figsize=(18,10))
    plots = []

    for i in range(6):
      plots.append(plt.subplot(2,3,i+1))

    plots[0].set_xlabel('p [hPa]')
    plots[1].set_xlabel('T [K]')
    plots[2].set_xlabel('RH')
 
    plots[3].set_xlabel('SO2  gas mole fraction')
    plots[4].set_xlabel('O3   gas mole fraction')  
    plots[5].set_xlabel('H2O2 gas mole fraction')

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
#      plots[3].plot(f.variables["SO2_g"][:]       , z, style[i])
#      plots[4].plot(f.variables["O3_g"][:]       , z, style[i])
#      plots[5].plot(f.variables["H2O2_g"][:] * 1000 , z, style[i])

    plt.savefig("plot_chem.svg")

