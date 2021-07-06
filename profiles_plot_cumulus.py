#import sys, os, subprocess
#sys.path.insert(0, "../../")
#sys.path.insert(0, "./")
#sys.path.insert(0, "/home/tomajev/singularity2/parcel/localinstall/lib/python3/dist-packages")

import sys, os, subprocess
sys.path.insert(0, "../")
sys.path.insert(0, "~/Piotr/IGF/parcel3/parcel/plots/comparison")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")
sys.path.insert(0, "/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")

from scipy.io import netcdf
import numpy as np
import math
from libcloudphxx import common

from parcel import parcel

def plot_profiles(fnc, output_folder="../outputs"):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties

    # ... plotting the results ...
    plt.figure(1, figsize=(18,10))
    plots    = []
    legend_l = []
    for i in range(6):
        plots.append(plt.subplot(2,3,i+1))

    plots[0].set_xlabel('p [hPa]')

    plots[1].ticklabel_format(useOffset=False) 
    plots[1].set_xlabel('th_d [K]')
    plots[2].set_xlabel('T [K]')
    plots[3].set_xlabel('kappa(rho_d :)) [kg/m3]')  
    plots[4].set_xlabel('rv [g/kg]')
    plots[5].set_xlabel('RH')

    for ax in plots:
        ax.set_ylabel('z [m]')

    z = fnc.variables["z"][:]
    plots[0].plot(fnc.variables["p"][:] / 100.   , z)
    plots[1].plot(fnc.variables["th_d"][:]       , z)
    plots[2].plot(fnc.variables["T"][:]          , z)
    plots[3].plot(fnc.variables["rhod"][:]       , z)
    plots[4].plot(fnc.variables["r_v"][:] * 1000 , z)
    plots[5].plot(
        fnc.variables["RH"][:]                     , z, 
        [fnc.variables["RH"][:].max()] * z.shape[0], z
        )
   
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "plot_profiles_onesim.svg"))

def main(dt=1):
    # running parcel model for different ways to solve for pressure  ...            
	RH_init = .98
	T_init = 280.
	p_init = 100000.
	r_init = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))
	outfile = "onesim_plot.nc"
	parcel(dt=1., outfreq = math.ceil(1./5),   outfile = outfile,\
                w = 5, T_0 = 280, p_0 = 100000, r_0 = r_init, z_max = 200, sstp_cond = 20, \
                sd_conc = 10000, \
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.03e-6, 0.14e-6], "gstdev": [1.28, 1.75], "n_tot": [90e6, 15e6]}}',
                out_bin = '{"cloud": {"rght": 2.5e-05, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-07}}'
                )
	fnc = netcdf.netcdf_file(outfile)
	plot_profiles(fnc)
	fnc.close()
#	subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()
