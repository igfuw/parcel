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

RH_init = .98
T_init  = 280.
p_init  = 100000.
r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))

w_list = [ 0.5, 1, 1.5, 2, 3, 4, 5]
sstp_cond = [1, 2, 5, 10, 20 ]

# initial values


for w in w_list:
    fig = plt.figure(figsize=(28,13))
    for st_cond in sstp_cond:
        outfile = "onesim_plot.nc"
        print "updraft velosity", w
        print "\nt condensation time step", st_cond
        outfile_nc = "timesteptest_cond=" + str(st_cond)+"updraft_velocity"+ str(w)+ ".nc"
        parcel(dt=1, outfreq = math.ceil(1./w),   outfile = outfile_nc,\
                w = w, T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 200, sstp_cond = st_cond, \
                sd_conc = 10000, \
                # aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.03e-6, 0.14e-6], "gstdev": [1.28, 1.75], "n_tot": [90e6, 15e6]}}',
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125e6, 65e6]}}',
                # out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0,1], "drwt": "wet", "nbin": 100, "lnli": "lin", "left": 1e-06 }}'
                out_bin = '{"cloud": {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 100, "moms": [0,1]}}')
                # '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125e6, 65e6]}}'
        output_folder="/home/piotr/Piotr/IGF/parcel/wyniki_momenty"
        fnc = netcdf.netcdf_file(outfile_nc)


        plots    = []
        legend_l = []
        z = fnc.variables["z"][:]
        num_cr = fnc.variables["number_of_rc"][:]
        sto = np.ones(len(z))
        ax = fig.add_subplot(141)
        ax.set_xlabel('RH [%]')
        ax.set_ylabel('Z [m]]')
        plt.grid()
        plt.plot(fnc.variables["RH"][:]*100 , z, label="1/"+str(st_cond))
        plt.plot(sto*100 , z,)
        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("Saturation for w= " + str(w)+ "[m/s]")
        ax = fig.add_subplot(142)
        ax.set_xlabel('Moment 0')
        ax.set_ylabel('Z [m]]')
        plt.grid()
        plt.plot(np.sum(fnc.variables["cloud_m0"][:, 1:],axis=1).tolist() , z, label="1/"+str(st_cond))
        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("0 moment for w= " + str(w) + "[m/s]")

        ax = fig.add_subplot(143)
        ax.set_xlabel('Critical radius')
        ax.set_ylabel('Z [m]]')
        plt.grid()
        plt.plot(num_cr  , z, label="1/"+str(st_cond))
        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("Critical radius for w= " + str(w) + "[m/s]")


        m0 = np.sum(fnc.variables["cloud_m0"][:],axis=1)
        m1 = np.sum(fnc.variables["cloud_m1"][:],axis=1)
        print m1/m0
        ax = fig.add_subplot(144)
        ax.set_xlabel('Mean radius')
        ax.set_ylabel('Z [m]]')
        plt.grid()
        plt.plot(m1/m0 , z, label="1/"+str(st_cond))
        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("Mean radius for w= " + str(w) + "[m/s]")
        if not os.path.exists(output_folder):
            subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "Condensation_time_step_and_Critical_radii_variation_for_cumulus2_updraft_velocity"+ str(w)+".svg"))

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
