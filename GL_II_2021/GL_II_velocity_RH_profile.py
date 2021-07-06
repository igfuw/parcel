import sys
sys.path.insert(0, "../")
# sys.path.insert(0, "/home/piotr/Piotr/IGF/parcel3/parcel/plots/comparison")
sys.path.insert(0, "/home/piotr-pc/Piotr/IGF/parcel3/parcel/plots/comparison")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")
# sys.path.insert(0, "/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
sys.path.insert(0, "/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")


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
# r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))
SD_conc = 500
H_max = 250
col = ['red', 'blue', 'green']

for j in range(len(mu)):
    fig = plt.figure(figsize=(25,18))
    i =0
    for w in w_list:

        print("updraft velosity", w, "mean value", mu[j])
        outfile_nc = "mean radius= "+str(mu[j])+",updraft_velocity= "+ str(w)+".nc"
        parcel(dt=1, outfreq = math.ceil(1./w),   outfile = outfile_nc,\
                # w = w, T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = H_max, sstp_cond = 10, \
                w = w, T_0 = T_init, p_0 = p_init, RH_0 = RH_init, z_max = H_max, sstp_cond = 10, \
                sd_conc = SD_conc, \

                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": ['+str(mu[j])+'], "gstdev": ['+str(stdev[j])+'], "n_tot": ['+str(N_tot[j])+']}}',
                # out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": '+str(mu[j]/10)+'}}'
                out_bin = '{"cloud": {"rght": 2.5e-04, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 1e-6}}'
                )
        output_folder="/home/piotr-pc/Piotr/IGF/parcel3/parcel/GL_II_2021/wyniki_momenty_TEST2"
        fnc = netcdf.netcdf_file(outfile_nc)


        plots    = []
        legend_l = []
        z = fnc.variables["z"][:]
        rhod = fnc.variables["rhod"][:]
        num_cr = fnc.variables["number_of_rc_m0"][:]
        mom0 = fnc.variables["cloud_m0"][:]
        # mom0 = np.sum(mom0,axis=1)*1.225
        mom3 = fnc.variables["cloud_m3"][:]
        # mom3 = np.sum(mom3,axis=1)
        sto = np.ones(len(z))
        mom0 = mom0.flatten()


        ax = fig.add_subplot(221)
        ax.set_xlabel('RH [%]')
        ax.set_ylabel('Z [m]]')
        plt.grid()
        plt.plot(fnc.variables["RH"][:]*100 , z, label=str(w), linewidth=2, color=col[i])
        plt.plot(sto*100 , z,)
        plt.legend(loc='upper left', title="updraft velocity [m/s]", frameon=1)
        plt.title("Relative humidity", loc="left", weight='bold')

        ax = fig.add_subplot(222)
        ax.set_xlabel(r"Concentration [$\frac{1}{cm^3}$]")
        ax.set_ylabel('Z [m]]')
        plt.grid()
        plt.plot( 1.225*mom0*1e-6 , z, label="cosnt rhod "+str(w), linewidth=2, color=col[i])
        plt.plot( rhod*mom0*1e-6 , z, label="rhod(H) "+str(w), linewidth=2, color=col[i])
        plt.legend(loc='upper left', title="updraft velocity [m/s]", frameon=1)
        plt.title("Concentration of " + r'$r>1\cdot 10^{-6} [m]$', loc="right", weight='bold')

        # ax = fig.add_subplot(223)
        # ax.set_xlabel('Concentration ' +r'$[\frac{1}{cm^{3}}]$')
        # ax.set_ylabel('Z [m]')
        # plt.grid()
        # plt.plot(num_cr*1e-6 , z, label=str(w))
        # plt.legend(loc='upper left', title="Cond sub-step", frameon=1)
        # plt.title("Concentration of " + r'$r>r_{crit.}$', weight='bold')

        ax = fig.add_subplot(223)
        ax.set_xlabel('LWC ' +r'$[\frac{g}{m^{3}}]$')
        ax.set_ylabel('Z [m]')
        plt.grid()
        plt.plot(mom3 * 4. / 3 * math.pi * 998.2 * 1e3 , z, label=str(w))
        # plt.plot(fnc.variables["cloud_m3"][:]*fnc.variables["cloud_m3"][:]*fnc.variables["cloud_m3"][:]*fnc.variables["cloud_m0"][:] * 4. / 3 * math.pi * 998.2 * 1e3 , z, label=str(w))
        plt.legend(loc='upper left', title="Cond sub-step", frameon=1)
        plt.title("Liquid Water Content", weight='bold')

        # m0 = np.sum(fnc.variables["cloud_m0"][:],axis=1)
        # m1 = np.sum(fnc.variables["cloud_m1"][:],axis=1)
        # m0 = fnc.variables["number_of_rc_m0"][:]
        # m1 = fnc.variables["number_of_rc_m1"][:]
        # ax = fig.add_subplot(224)
        # ax.set_xlabel('mean radius ' + r'$[\mu m]$')
        # ax.set_ylabel('Z [m]')
        # plt.grid()
        # plt.ylim((0, 260))
        # plt.plot(m1/m0*1e6, z, label=str(w))
        # plt.legend(loc='upper left', title="Cond sub-step", frameon=1)
        # plt.title("Mean critical radius ", weight='bold')

        ax = fig.add_subplot(224)
        ax.set_xlabel('mean volume radius ' + r'$[\mu m]$')
        ax.set_ylabel('Z [m]')
        plt.grid()
        plt.ylim((0, 260))
        plt.plot(mom3*1e6, z, label=str(w))
        plt.legend(loc='upper left', title="Cond sub-step", frameon=1)
        plt.title("Mean volume radius ", weight='bold')

        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        plt.suptitle("$\mu$= "+str(mu[j]*1e6) + "[μm], $\sigma$ = " + str(stdev[j]) +", N = "+str(N_tot[j]/1e6) +r"[$\frac{\#\cdot 1e6}{m^{-3}}$], num of SD = "+str(SD_conc), y=0.995, weight='bold')
        plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=1.0)
        i +=1
        if not os.path.exists(output_folder):
            subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "mean radius= "+str(mu[j])+"Number of Super_droplets= "+str(SD_conc)+ " LWC.png"))
