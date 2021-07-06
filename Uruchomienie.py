import sys
sys.path.insert(0, "../")
sys.path.insert(0, "/home/piotr/Piotr/IGF/parcel/plots/comparison")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")

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

w_list = [0.1, 0.5, 1, 2, 5]
sstp_cond = [1, 2, 5, 10, 20]

# initial values
RH_init = .98
T_init  = 280.
p_init  = 100000.
r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))

for w in w_list:
    fig = plt.figure(figsize=(28,13))
    for st_cond in sstp_cond:
        outfile = "onesim_plot.nc"
        print("updraft velosity", w)
        print("\nt condensation time step", st_cond)
        outfile_nc = "timesteptest_cond=" + str(st_cond)+"updraft_velocity"+ str(w)+ ".nc"
        parcel(dt=1, outfreq = math.ceil(1./w),   outfile = outfile_nc,\
                w = w, T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 200, sstp_cond = st_cond, \
                sd_conc = 10000, \
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.03e-6, 0.14e-6], "gstdev": [1.28, 1.75], "n_tot": [180e6, 30e6]}}',
                out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-07}}'
                )
        output_folder="/home/piotr/Piotr/IGF/parcel3/parcel/wyniki_momenty_R2"
        fnc = netcdf.netcdf_file(outfile_nc)
