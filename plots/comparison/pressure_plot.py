import sys, os, subprocess
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "../../")
from scipy.io import netcdf
import numpy as np

from libcloudphxx import common
from parcel import parcel

pprof_list = ["pprof_const_rhod", "pprof_const_th_rv", "pprof_piecewise_const_rhod"]

def plot_pressure_opt(fnc, output_folder="../outputs"):
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
    # the different ways of solving for pressure come from different assumptions about the density profile
    # but those assumptions are just used when calculating the profile of pressure
    # later on the rho_d profile can be calculated (and is not the same as the one assumed)
    # so the kappa here is the actual profile of rho_d during the simulation (different than the one assumed)
    plots[3].set_xlabel('kappa(rho_d :)) [kg/m3]')  
    plots[4].set_xlabel('rv [g/kg]')
    plots[5].set_xlabel('RH')

    for ax in plots:
        ax.set_ylabel('z [m]')

    style = ["g.-", "b.-","r.-"]
    for i, pprof in enumerate(pprof_list):
        z = fnc[pprof].variables["z"][:]
        plots[0].plot(fnc[pprof].variables["p"][:] / 100.   , z, style[i])
        plots[1].plot(fnc[pprof].variables["th_d"][:]       , z, style[i])
        plots[2].plot(fnc[pprof].variables["T"][:]          , z, style[i])
        plots[3].plot(fnc[pprof].variables["rhod"][:]       , z, style[i])
        plots[4].plot(fnc[pprof].variables["r_v"][:] * 1000 , z, style[i])
        plots[5].plot(
              fnc[pprof].variables["RH"][:]                     , z, style[i], 
                [fnc[pprof].variables["RH"][:].max()] * z.shape[0], z, style[i]
        )
        legend_l.append(pprof)
    plots[0].legend(legend_l, loc=1, prop = FontProperties(size=10))
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "plot_pressure.svg"))

def main(dt=0.1):
    # running parcel model for different ways to solve for pressure  ...            
    fnc = {}
    for pprof in pprof_list:
        outfile = "test_" + pprof + ".nc"
        parcel(dt=dt, outfreq = 100, pprof = pprof, outfile=outfile)
        fnc[pprof] = netcdf.netcdf_file(outfile)
    plot_pressure_opt(fnc)
    subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()
