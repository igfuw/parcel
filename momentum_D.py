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
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125e6, 65e6]}}',
                out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-07}}'
                )
        output_folder="/home/piotr/Piotr/IGF/parcel3/parcel/wyniki_momenty_D"
        fnc = netcdf.netcdf_file(outfile_nc)


        plots    = []
        legend_l = []
        z = fnc.variables["z"][:]
        sto = np.ones(len(z))
        ax = fig.add_subplot(121)
        ax.set_xlabel('RH [%]')
        ax.set_ylabel('Z [m]]')

        plt.grid()
        plt.plot(fnc.variables["RH"][:]*100 , z, label="1/"+str(st_cond))
        plt.plot(sto*100 , z,)

        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("Saturation for updraft velocity= " + str(w)+ "[m/s]")
        ax = fig.add_subplot(122)
        ax.set_xlabel('Moment 0')
        ax.set_ylabel('Z [m]]')
        plt.grid()
        plt.plot(fnc.variables["cloud_m0"][:] , z, label="1/"+str(st_cond))
        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("0 moment for updraft velocity= " + str(w) + "[m/s]")
        if not os.path.exists(output_folder):
            subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "Condensation_time_step_variation_for_stratocumulus_updraft_velocity"+ str(w)+".svg"))

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
