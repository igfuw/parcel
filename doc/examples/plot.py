import sys, os, subprocess
sys.path.insert(0, "../../")
sys.path.insert(0, "./")

#listing00
from scipy.io import netcdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from parcel import parcel

# define parcel arguments
outfile = "example_py.nc"
dt      = 0.5
outfreq = 10

#run parcel model
parcel(dt=dt, outfreq = outfreq, outfile=outfile)

# open ncdf file with model results
ncfile = netcdf.netcdf_file(outfile)

# plot the results
plt.figure(1, figsize=(20,10))
plots    = []
legend_l = []
for i in range(3):
    plots.append(plt.subplot(1,3,i+1))

plots[0].set_xlabel('p [hPa]')
plots[1].set_xlabel('T [K]')
plots[2].set_xlabel('RH')

for ax in plots:
    ax.set_ylabel('z [m]')

z = ncfile.variables["z"][:]
plots[0].plot(ncfile.variables["p"][:] / 100. , z)
plots[1].plot(ncfile.variables["T"][:]        , z)
plots[2].plot(ncfile.variables["RH"][:]       , z) 
   
plt.savefig("doc_python.pdf")
#listing01

