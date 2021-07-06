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

w_list = [1, 3, 5]
sstp_cond = [10]
mu = [0.24e-6, 0.06e-6, 0.54e-6, 0.14e-6]
N_tot = [90e6, 70e6, 90e6, 80e6]
stdev = [1.58, 1.38, 2.2, 1.2]

# initial values
RH_init = .98
T_init  = 288.
p_init  = 100000.
r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))
SD_conc = 10000
col = ['red', 'blue', 'yellow']

for j in range(len(mu)):
    fig = plt.figure(figsize=(28,13))
    i =0
    for w in w_list:

        outfile = "onesim_plot.nc"
        print("updraft velosity", w, "mean value", mu[j])
        outfile_nc = "TEST_04_05_cond=10 updraft_velocity"+ str(w)+"mean radius"+str(mu[j]) ".nc"
        parcel(dt=1, outfreq = math.ceil(1./w),   outfile = outfile_nc,\
                w = w, T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 200, sstp_cond = 10, \
                sd_conc = SD_conc, \

                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": ['+str(mu[j])+'], "gstdev": ['+str(stedv[j])+'], "n_tot": ['+str(N_tot[j])+']}}',
                out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 1e-6}}'
                )
    output_folder="/home/piotr/Piotr/IGF/parcel3/parcel/GL_II_2021/wyniki_momenty_TEST"
    fnc = netcdf.netcdf_file(outfile_nc)


    plots    = []
    legend_l = []
    z = fnc.variables["z"][:]
    sto = np.ones(len(z))
    ax = fig.add_subplot(122)
    ax.set_xlabel('RH [%]')
    ax.set_ylabel('Z [m]]')

    plt.grid()
    plt.plot(fnc.variables["RH"][:]*100 , z, label=str(w), linewidth=2, color=col[i])
    plt.plot(sto*100 , z,)

    plt.legend(loc='upper left', title="updraft velocity [m/s]", frameon=1)
    plt.title("Condensation time step 20 ")
    ax = fig.add_subplot(121)
    ax.set_xlabel('Concentration [cm^-3]')
    ax.set_ylabel('Z [m]]')
    plt.grid()
    plt.plot(fnc.variables["cloud_m0"][:]*1e-6 , z, label=str(w), linewidth=2, color=col[i])
    plt.legend(loc='upper left', title="updraft velocity [m/s]", frameon=1)
    plt.title("Condensation time step 10")
    i +=1
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
plt.savefig(os.path.join(output_folder, "TEST_04_05, 0.5um, sstp=10, sd_conc= "+str(SD_conc)+".png"))
