import sys, os, subprocess
sys.path.insert(0, "../../")
sys.path.insert(0, "./")

from scipy.io import netcdf
import numpy as np

from parcel import parcel

def thesis_profiles(fnc, output_folder="../outputs"):

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties

    # ... plotting the results ...
    plt.figure(1)
    plt.figure(figsize=(28,13))
    plt.rcParams.update({'font.size': 30})
    plots    = []
    legend_l = []
    for i in range(3):
        plots.append(plt.subplot(1,3,i+1))

    for ax in plots:
      ax.set_ylabel('z [m]')
      ax.set_ylim([0, 1400])
      ax.set_yticks([0, 200, 400, 600, 800, 1000, 1200, 1400])
      ax.tick_params(axis='x', pad=15)
      ax.tick_params(axis='y', pad=15)

    plots[0].ticklabel_format(useOffset=False) 
    plots[0].set_xlabel('$\mathrm{\Theta_d}$ [K]')
    plots[1].set_xlabel('$\mathrm{r_v}$ [g/kg]')
    plots[2].set_xlabel('RH [%]')

    z = fnc.variables["z"][:]
    plots[0].plot( fnc.variables["th_d"][:]                  , z, "b.-", ms=15, lw=4.)
    plots[1].plot( fnc.variables["r_v"][:] * 1000            , z, "b.-", ms=15, lw=4.)
    plots[2].plot( fnc.variables["RH"][:]                    , z, "b.-", ms=15, lw=4.)
    plots[2].plot([fnc.variables["RH"][:].max()] * z.shape[0], z, "r.-", ms=15, lw=4.)
    plots[2].plot([1] * z.shape[0]                           , z, "g.-", ms=15, lw=4.)
   
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "thesis_profiles.svg"))

def main(dt=1):
    # running parcel model for different ways to solve for pressure  ...            
    outfile = "thesis_plot.nc"
    parcel(dt=dt, outfreq = 10, outfile=outfile)
    fnc = netcdf.netcdf_file(outfile)
    thesis_profiles(fnc)
    fnc.close()
    subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()
