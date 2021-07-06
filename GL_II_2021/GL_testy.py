import sys
sys.path.insert(0, "../")
sys.path.insert(0, "/home/piotr/Piotr/IGF/parcel3/parcel/plots/comparison")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")
sys.path.insert(0, "/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")


from parcel import parcel
from libcloudphxx import common
from timestep_plot import timestep_plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.io import netcdf
import numpy as np
import pytest
import os, glob
import subprocess
import pdb
import math
plt.rcParams.update({'font.size': 18})

outfile_nc ="mean radius= 2.4e-07,updraft_velocity= 1.nc"
fnc = netcdf.netcdf_file(outfile_nc)

plots    = []
legend_l = []
z = fnc.variables["z"][:]
num_cr = fnc.variables["number_of_rc_m0"][:]
mom0 = fnc.variables["cloud_m0"][:]
mom0 = np.sum(mom0,axis=1)
sto = np.ones(len(z))

print(4/3)
# plt.plot(mom0, z)
# plt.savefig("LWC.png")
