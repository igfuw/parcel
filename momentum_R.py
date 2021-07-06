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
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.animation import FuncAnimation
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
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
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.03e-6, 0.14e-6], "gstdev": [1.28, 1.75], "n_tot": [90e6, 15e6]}}',
                out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-07}}'
                )
        output_folder="/home/piotr/Piotr/IGF/parcel3/parcel/wyniki_momenty_R"
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
    plt.savefig(os.path.join(output_folder, "Condensation_time_step_variation_for_cumulus_updraft_velocity"+ str(w)+".svg"))

        # print np.percentile(fnc.variables["cloud_r_wet"][:], 10)*1e6
        # print np.percentile(fnc.variables["cloud_r_wet"][:], 25)*1e6
        # print np.percentile(fnc.variables["cloud_r_wet"][:], 50)*1e6
        # print np.percentile(fnc.variables["cloud_r_wet"][:], 75)*1e6
        # print np.percentile(import matplotlib.patches as patches

    '''
        # print np.percentile(fnc.variables["cloud_m0"][0],10)
        ax = fig.add_subplot(131)
        ax.set_xlabel('RH [%]')
        ax.set_ylabel('Z [m]')

        plt.grid()
        plt.plot(fnc.variables["RH"][:]*100 , z, label="1/"+str(st_cond))
        plt.plot(sto*100 , z,)

        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("Saturation for updraft velocity= " + str(w)+ "[m/s]")
        ax = fig.add_subplot(132)
        # ax.set_xlabel('Moment 0')
        # ax.set_ylabel('Z [m]]')
        ax.set_xlabel('Radius')
        ax.set_ylabel('number of particels ')
        plt.grid()
        # plt.plot(fnc.variables["cloud_r_wet"][:]*1e6, fnc.variables["cloud_m0"][1])
        plt.plot(fnc.variables["cloud_m0"][:].sum(axis=1) , z, label="1/"+str(st_cond))
        plt.legend(loc='upper left', title="Condensation time step", frameon=1)
        plt.title("0 moment for updraft velocity= " + str(w) + "[m/s]")



        # ax = fig.add_subplot(133)
        # n, bins, patches = plt.hist(fnc.variables["cloud_m0"][55], 10, normed=1, facecolor='green', alpha=0.75)
        # plt.hist(fnc.variables["cloud_m0"][10], bins=26, label="100 1/"+str(st_cond))
        # plt.legend(loc='upper right', title="Condensation time step", frameon=1)
        # plt.hist(fnc.variables["cloud_m0"][-1], bins=26)




        ax = fig.add_subplot(133)
        # ax.set_xlabel('Radius')
        # ax.set_ylabel('Number')
        # plt.grid()
        # axs[0].hist(x, bins=n_bins)
        # plt.hist(fnc.variables["cloud_m0"][-1], bins=26)
        # plt.hist(fnc.variables["cloud_m0"][1], bins=26)
        # suma = np.ones(len(z))
        # for i in range (len(z)):
            # suma[i] = fnc.variables["cloud_m0"][i].sum(axis=0)
        # suma1 = fnc.variables["cloud_m0"][100].sum(axis=0)
        # suma2 = fnc.variables["cloud_m0"][-1].sum(axis=0)

        # print suma[-1]

        # print np.cumsum(fnc.variables["cloud_m0"][100]/suma1)
        # print np.cumsum(fnc.variables["cloud_m0"][-1]/suma2)
        # plt.plot(fnc.variables["cloud_r_wet"][:]*1e6, np.cumsum(fnc.variables["cloud_m0"][100]/suma1),label="100 1/"+str(st_cond))
        # plt.plot(fnc.variables["cloud_r_wet"][:]*1e6, np.cumsum(fnc.variables["cloud_m0"][-1]/suma2), label="110 1/"+str(st_cond))
        # plt.plot( np.percentile(fnc.variables["cloud_m0"][1], 25)*1e6, c='y')
        # plt.plot( np.percentile(fnc.variables["cloud_m0"][1], 50)*1e6, c='b')
        # plt.plot( np.percentile(fnc.variables["cloud_m0"][1], 75)*1e6, c='g')
        plt.plot( fnc.variables["r_v"][:], z,  c='m')
        # plt.legend(loc='lower right', title="Condensation time step", frameon=1)
        # plt.title("0 moment for updraft velocity= " + str(w) + "[m/s]")
        # fig = plt.figure()
        # ax = plt.axes(xlim=(0, 25), ylim=(0, 1))
        # line, = ax.plot([], [], lw=3)
        #
        # def init():
        #     line.set_data([], [])
        #     return line,
        # def animate(i):
        #     x = fnc.variables["cloud_r_wet"][:]*1e6
        #     y = np.cumsum(fnc.variables["cloud_m0"][i]/fnc.variables["cloud_m0"][i].sum(axis=0))
        #     line.set_data(x, y)
        #     return line,
        #
        # anim = FuncAnimation(fig, animate, init_func=init,
        #                        frames=200, interval=1, blit=True)

        if not os.path.exists(output_folder):
            subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "Condensation_time_step_variation_for_cumulus_updraft_velocity"+ str(w)+".svg"))
    # anim.save('sine_wave.html', writer='imagemagick')
    #
    # data = {"sstp_cond" : sstp_cond, "w" : w_list}
    '''

# ... plotting the results ...


#         RH_max = f_out.RH_max
#         N_end  = f_out.variables["radii_m0"][-1,0] # concentration of drops > 1e-6 m
#
#         RH_list.append((RH_max - 1)*100)  # [%]
#         N_list.append(N_end / 1e6)        # [1/mg]
#
# data = {"RH" : RH_list, "N" : N_list, "dt" : Dt_list}
