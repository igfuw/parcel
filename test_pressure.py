from libcloudphxx import common
import matplotlib.pyplot as plt
from scipy.io import netcdf
from parcel import parcel, Pprof

# big timestep to see some differences
dt = 5

# running parcel model for different ways to solve for pressure  ...
parcel(dt=dt, pprof=Pprof.hydro_const_rhod,            outfile="test_hydro_const_rhod.nc")
parcel(dt=dt, pprof=Pprof.hydro_const_th_rv,           outfile="test_hydro_const_th_rv.nc")
parcel(dt=dt, pprof=Pprof.hydro_piecewise_const_th_rv, outfile="test_hydro_piecewise_const_th_rv.nc")
parcel(dt=dt, pprof=Pprof.hydro_old_drops,             outfile="test_hydro_old_drops.nc")

# ... and the rest is just for plotting
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
