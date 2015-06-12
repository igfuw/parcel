import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
import matplotlib.pyplot as plt
from scipy.io import netcdf
from parcel import parcel, Pprof
import numpy as np
import pytest

@pytest.mark.parametrize("dt", [1])
def test_pressure(dt):
    # running parcel model for different ways to solve for pressure  ...
    parcel(dt=dt, pprof=Pprof.hydro_const_rhod,            outfile="test_hydro_const_rhod.nc")
    parcel(dt=dt, pprof=Pprof.hydro_const_th_rv,           outfile="test_hydro_const_th_rv.nc")
    parcel(dt=dt, pprof=Pprof.hydro_piecewise_const_th_rv, outfile="test_hydro_piecewise_const_th_rv.nc")
    parcel(dt=dt, pprof=Pprof.hydro_old_drops,             outfile="test_hydro_old_drops.nc")

    # ... plotting the results ...
    f_out = {
      "WWG-LPW"   : netcdf.netcdf_file("test_hydro_const_rhod.nc", "r"),
      "icicle"    : netcdf.netcdf_file("test_hydro_const_th_rv.nc", "r"),
      "piecewise" : netcdf.netcdf_file("test_hydro_piecewise_const_th_rv.nc", "r"),
      "drops"     : netcdf.netcdf_file("test_hydro_old_drops.nc", "r")
    }

    style = {
      "WWG-LPW"   : "g.-",
      "icicle"    : "b.-",
      "piecewise" : "r.-",
      "drops"     : "m.-"
    }

    plt.figure(1, figsize=(18,10))
    plots = []

    for i in range(6):
      plots.append(plt.subplot(2,3,i+1))

    plots[0].set_xlabel('p [hPa]')

    plots[1].ticklabel_format(useOffset=False) 
    plots[1].set_xlabel('th_d [K]')
    plots[2].set_xlabel('T [K]')
    # the different ways of solving for pressure come from different assumptions about the density profile
    # but those assumptions are just used when calculating the profile of pressure
    # later on the rho_d profile can be calculated (and is not the same as the one assumed)
    # so the kappa here is the actual profile of rho_d during the simulation (different than the one assumed)
    plots[3].set_xlabel('kappa(rho_d :)) [kg/m3]')  
    plots[4].set_xlabel('rv [g/kg]')
    plots[5].set_xlabel('RH')

    for ax in plots:
      ax.set_ylabel('z [m]')

    for i, f in f_out.iteritems():
      z = f.variables["z"][:]
      plots[0].plot(f.variables["p"][:] / 100.   , z, style[i])
      plots[1].plot(f.variables["th_d"][:]       , z, style[i])
      plots[2].plot(f.variables["T"][:]          , z, style[i])
      plots[3].plot(f.variables["rhod"][:]       , z, style[i])
      plots[4].plot(f.variables["r_v"][:] * 1000 , z, style[i])
      plots[5].plot(
	f.variables["RH"][:]                     , z, style[i], 
	[f.variables["RH"][:].max()] * z.shape[0], z, style[i]
      )

    plt.savefig("plot.svg")

    # ... and checking wheather those four different methods dont differ too much
    z_all = p_all = thd_all = T_all = rhod_all = rv_all = RH_max = [] 
    for i, f in f_out.iteritems():
      z_all    = np.append(z_all,    f.variables["z"][-1])
      p_all    = np.append(p_all,    f.variables["p"][-1])
      thd_all  = np.append(thd_all,  f.variables["th_d"][-1])
      T_all    = np.append(T_all,    f.variables["T"][-1])
      rhod_all = np.append(rhod_all, f.variables["rhod"][-1])
      rv_all   = np.append(rv_all,   f.variables["r_v"][-1])
      RH_max   = np.append(RH_max,   f.variables["RH"][:].max()) 

    assert (z_all == 200).all()
    assert (p_all >= 99039).all     and (p_all <= 99338).all
    assert (thd_all >= 302.50).all  and (thd_all <= 302.66).all
    assert (T_all >= 299.14).all    and (T_all <= 299.23).all
    assert (rhod_all >= 1.118).all  and (rhod_all <= 1.121).all
    assert (rv_all >= 0.02168).all  and (rv_all <= 0.02174).all
    assert (RH_max >= 1.0040).all   and (RH_max <= 1.0047).all
