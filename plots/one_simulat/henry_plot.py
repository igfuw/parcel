import sys
import os
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "../../")

from scipy.io import netcdf
import numpy as np
import pytest
import math
import subprocess

from parcel import parcel
from libcloudphxx import common as cm
import functions as fn

def plot_henry(data, output_folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    chem = ["SO2", "O3", "H2O2", "CO2", "HNO3", "NH3"]

    z    = data.variables["z"][:]
    t    = data.variables["t"][:]
    vol  = np.squeeze(data.variables["radii_m3"][:]) * 4/3. * math.pi
    T    = data.variables["T"][:]
    p    = data.variables["p"][:]
    rhod = data.variables["rhod"][:]
    conc_H = np.squeeze(data.variables["chem_H"][:]) / cm.M_H / vol

    plt.figure(1, figsize=(18,14))
    plots    = []
    legend_l = []

    for i in range(6):
        plots.append(plt.subplot(2,3,i+1))
        lab = chem[i] + " in droplets [kg/kg dry air]"
        plots[i].set_xlabel(lab)

    for ax in plots:
        ax.set_ylabel('t [s]')

    for i in range(6):
        mixr_g = data.variables[chem[i]+"_g"][:]
        plots[i].plot(fn.henry_teor(chem[i], p, T, vol, mixr_g, rhod, conc_H), t, "r.-", label="Henry(T, pH)")
        plots[i].plot(data.variables[chem[i]+"_a"][:], t, "b.-", label="in drop")
        plots[i].legend(loc='upper left')
        plots[i].ticklabel_format(axis = 'x', style='sci', scilimits=(-2,2))

    plt.savefig(os.path.join(output_folder, "plot_Henry.pdf"))
    plt.clf()

def main():
    RH_init = .99999
    T_init  = 300.
    p_init  = 100000.
    r_init  = cm.eps * RH_init * cm.p_vs(T_init) / (p_init - RH_init * cm.p_vs(T_init))

    # calculate rhod for initial gas mixing ratio
    th_0      = T_init * (cm.p_1000 / p_init)**(cm.R_d / cm.c_pd)
    rhod_init = cm.rhod(p_init, th_0, r_init)

    SO2_g_init  = fn.mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
    O3_g_init   = fn.mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
    H2O2_g_init = fn.mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
    CO2_g_init  = fn.mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
    NH3_g_init  = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
    HNO3_g_init = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

    outfreq     = 50
    z_max       = 50.
    outfile     = "test_chem_henry.nc"
    dt          = 0.1
    wait        = 1000
    sd_conc     = 1

    # run parcel run!
    parcel(dt = dt, z_max = z_max, outfreq = outfreq, wait=wait,\
            T_0 = T_init, p_0 = p_init, r_0 = r_init, sd_conc=sd_conc,\
            SO2_g = SO2_g_init, O3_g  = O3_g_init,  H2O2_g = H2O2_g_init,\
            CO2_g = CO2_g_init, NH3_g = NH3_g_init, HNO3_g = HNO3_g_init,\
            outfile = outfile,\
            chem_dsl = True, chem_dsc = False, chem_rct = False,\
            out_bin = \
            '{"radii": {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 3]},\
              "chem" : {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1,\
                  "moms": ["O3_a", "H2O2_a", "SO2_a", "CO2_a", "NH3_a", "HNO3_a", "H"]}}')

    data = netcdf.netcdf_file(outfile,"r")

    #plot
    plot_henry(data, output_folder = "../outputs")

    #cleanup
    data.close()
    subprocess.call(["rm", "test_chem_henry.nc"])

if __name__ == '__main__':
    main()
