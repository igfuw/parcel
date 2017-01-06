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
    plt.figure(1, figsize=(12,8))
    plots    = []
    legend_l = []
    for i in range(3):
        plots.append(plt.subplot(1,3,i+1))

    plots[0].ticklabel_format(useOffset=False) 
    plots[0].set_xlabel('th_d [K]')
    plots[1].set_xlabel('rv [g/kg]')
    plots[2].set_xlabel('RH')

    for ax in plots:
        ax.set_ylabel('z [m]')

    z = fnc.variables["z"][:]
    plots[0].plot(fnc.variables["th_d"][:]       , z)
    plots[1].plot(fnc.variables["r_v"][:] * 1000 , z)
    plots[2].plot(
        fnc.variables["RH"][:]                     , z, 
        [fnc.variables["RH"][:].max()] * z.shape[0], z,
        [1] * z.shape[0],                            z
        )
   
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
