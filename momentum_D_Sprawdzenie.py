import sys
sys.path.insert(0, "../")
sys.path.insert(0, "~/Piotr/IGF/parcel3/parcel/plots/comparison")
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

w_list = [4, 4.5, 5]
# w_list = [10]
# sstp_cond = [1, 2, 5, 10, 20]
#SPRAWDZENIE
# sstp_cond = [20]

# initial values
RH_init = .98
T_init  = 280.
p_init  = 100000.
r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))
SD_conc = 10000
col = ['red', 'blue', 'yellow']

fig = plt.figure(figsize=(28,13))
i =0
for w in w_list:

    outfile = "onesim_plot.nc"
    print("updraft velosity", w)
    # print("\nt condensation time step", st_cond)
    outfile_nc = "TEST_timesteptest_cond=20 updraft_velocity"+ str(w)+ ".nc"
    parcel(dt=1, outfreq = math.ceil(1./w),   outfile = outfile_nc,\
                w = w, T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 200, sstp_cond = 20, \
                sd_conc = SD_conc, \
                #DYCOMS
                # aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125e6, 65e6]}}',
                #RICO
                # aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.03e-6, 0.14e-6], "gstdev": [1.28, 1.75], "n_tot": [90e6, 15e6]}}',
                #RICOx2
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.03e-6, 0.14e-6], "gstdev": [1.28, 1.75], "n_tot": [180e6, 30e6]}}',
                out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-6}}'
                )
                #5e-07
    output_folder="/home/piotr/Piotr/IGF/parcel3/parcel/wyniki_momenty_D"
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
    plt.title("Condensation time step 20")
    i +=1
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
plt.savefig(os.path.join(output_folder, "Ricox2_test2, 0.5um, sstp=20, sd_conc= "+str(SD_conc)+".png"))

    #     f_out  = netcdf.netcdf_file(outfile_nc, "r")
    # data = {"sstp_cond" : sstp_cond, "w" : w_list}


# ... plotting the results ...


#         RH_max = f_out.RH_max
#         N_end  = f_out.variables["radii_m0"][-1,0] # concentration of drops > 1e-6 m
#
#         RH_list.append((RH_max - 1)*100)  # [%]
#         N_list.append(N_end / 1e6)        # [1/mg]
#
# data = {"RH" : RH_list, "N" : N_list, "dt" : Dt_list}
