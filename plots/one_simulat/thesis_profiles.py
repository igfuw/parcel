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
    #from mpl_toolkits.axes_grid.inset_locator import inset_axes

    z = fnc.variables["z"][:]

    # ... plotting the results ...
    fig = plt.figure(figsize=(28,13))
    plt.rcParams.update({'font.size': 30})
  
    ax = fig.add_subplot(131) 
    ax.set_ylim([0, 1300])
    ax.set_yticks([  0,100, 200, 300, 400, 500, 600, 700,    800,  900, 1000, 1100, 1200, 1300])
    y_labels=[    -200,  0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400]
    ax.set_yticklabels(y_labels)
    ax.tick_params(axis='x', pad=15)
    ax.tick_params(axis='y', pad=15)
    #ax.ticklabel_format(useOffset=False) 
    ax.set_xlabel('$\mathrm{\Theta_d}$ [K]')
    ax.set_xlim([289, 297])
    ax.set_xticks([289, 291, 293, 295, 297]) 
    ax.set_ylabel('time above cloud base [s]')
    plt.grid()
    plt.plot( fnc.variables["th_d"][:]                  , z, "b.-", ms=15, lw=4.)

    ax = fig.add_subplot(132) 
    ax.set_ylim([0, 1300])
    ax.set_yticks([0,100, 200,300, 400,500, 600,700, 800,900, 1000,1100, 1200, 1300])
    ax.tick_params(axis='x', pad=15)
    ax.tick_params(axis='y', pad=15)
    ax.set_xlabel('$\mathrm{r_v}$ [g/kg]')
    plt.grid()
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    plt.plot( fnc.variables["r_v"][:] * 1000            , z, "b.-", ms=15, lw=4.)

    ax = fig.add_subplot(133) 
    ax.set_ylim([0, 1300])
    ax.set_yticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300])
    y_labels=[     0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,    1,  1.1,  1.2,  1.3]
    ax.set_yticklabels(y_labels)
    ax.yaxis.tick_right()
    ax.set_ylabel('height above the ground [km]')
    ax.yaxis.set_label_position("right")
    ax.tick_params(axis='x', pad=15)
    ax.tick_params(axis='y', pad=15)
    ax.set_xlabel('RH')
    ax.set_xlim([0.95, 1.01])
    ax.set_xticks([0.95, 0.97, 0.99, 1.01]) 
    plt.grid()
    plt.plot( fnc.variables["RH"][:]                    , z, "b.-", ms=15, lw=4.)
    #plt.plot([fnc.variables["RH"][:].max()] * z.shape[0], z, "r.-", ms=15, lw=4.)
    plt.plot([1] * z.shape[0]                           , z, "r-", lw=4.)

    #inset_axes = inset_axes(ax, 
    #                width="50%", # width = 30% of parent_bbox
    #                height=4.0, # height : 1 inch
    #                loc=2)
    #plt.plot( fnc.variables["RH"][2:6]                         , z[2:6], "b.-", ms=15, lw=4.)
    #plt.plot([fnc.variables["RH"][2:6].max()] * z[2:6].shape[0], z[2:6], "r.-", ms=15, lw=4.)
    #plt.plot([1] * z[2:6].shape[0]                             , z[2:6], "g.-", ms=15, lw=4.)
   
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "Kreidenweis_thesis_profiles.pdf"))

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
