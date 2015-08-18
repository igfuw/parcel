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

# TODO - don't repeat with the same function in unit_test/test_chem_henry
def henry_teor(chem, p, T, vol, conc_g):
    H   = getattr(cm, "H_"  +chem)
    dHR = getattr(cm, "dHR_"+chem)
    if chem in ["SO2", "CO2", "NH3"]:
        molar_mass = getattr(cm, "M_"+chem+"_H2O")
    else:
        molar_mass = getattr(cm, "M_"+chem)

    # correction to Henry const. due to temperature
    henry_T = H * np.exp(-1 * dHR * (1./T - 1./298))

    # dissolved  = partial prsessure * Henry_const * molar mass * drop volume
    henry_exp    = conc_g * p * henry_T * molar_mass * vol
    return henry_exp

def plot_henry(data, chem_sys, output_folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    chem = ["SO2", "O3", "H2O2", "CO2", "HNO3", "NH3"]

    z    = data.variables["z"][:]
    t    = data.variables["t"][:]
    vol  = np.squeeze(data.variables["radii_m3"][:]) * 4/3. * math.pi
    lwc  = vol * 998.2 * 1000.
    T    = data.variables["T"][:]
    p    = data.variables["p"][:]

    plt.figure(1, figsize=(18,14))
    plots    = []
    legend_l = []

    for i in range(6):
        plots.append(plt.subplot(4,3,i+1))
        lab = chem[i] + "dissolved in droplets kg/kg dry air"
        plots[i].set_xlabel(lab)

    plots.append(plt.subplot(4,3,7))
    plots[6].set_xlabel('z [m]')
    plots.append(plt.subplot(4,3,8))
    plots[7].set_xlabel('p [hPa]')
    plots.append(plt.subplot(4,3,9))
    plots[8].set_xlabel('rv [g/kg]')
    plots.append(plt.subplot(4,3,10))
    plots[9].set_xlabel('gas vol.conc SO2 [ppb]')
    plots[9].set_xticks([0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    plots[9].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3'])
    plots[9].set_xlim([0., 0.3])
    plots.append(plt.subplot(4,3,11))
    plots[10].set_xlabel('lwc [g/kg]')
    plots.append(plt.subplot(4,3,12))
    plots[11].set_xlabel('T [K]')

    for ax in plots:
        ax.set_ylabel('t [s]')

    for i in range(6):
        conc_g = data.variables[chem[i]+"_g"][-1]
        plots[i].plot(henry_teor(chem[i], p, T, vol, conc_g), t, "r.-", label="Henry")
        plots[i].plot(data.variables[chem[i]+"_a"][:], t, "b.-", label="in drop")
        plots[i].legend(loc='upper left')

    plots[6].plot(z, t)
    plots[7].plot(p * 1000, t)
    plots[8].plot(data.variables["r_v"][:] * 1000, t)
    plots[9].plot(data.variables["SO2_g"][:] * 1e9, t)
    plots[10].plot(lwc, t)
    plots[11].plot(T, t)

    plt.savefig(os.path.join(output_folder, "plot_Henry_"+chem_sys+".pdf"))
    plt.clf()

def main():
    RH_init = .99999
    T_init  = 300.
    p_init  = 100000.
    r_init  = cm.eps * RH_init * cm.p_vs(T_init) / (p_init - RH_init * cm.p_vs(T_init))

    SO2_g_init  = 200e-12
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12
    outfreq     = 50
    z_max       = 50.
    outfile     = "test_chem_dsl.nc"
    dt          = 0.1
    wait        = 1000

    # run parcel model for open chem system  and only for dissolving chem species into droplets
    parcel(dt = dt, z_max = z_max, outfreq = outfreq, wait=wait,\
            T_0 = T_init, p_0 = p_init, r_0 = r_init,\
            SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
            CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
            chem_sys = 'open',   outfile = outfile,\
            chem_dsl = True, chem_dsc = False, chem_rct = False,\
            out_bin = \
            '{"radii": {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 3]},\
              "chem" : {"rght": 1, "left": 0, "drwt": "wet", "lnli": "lin", "nbin": 1,\
                  "moms": ["O3_a", "H2O2_a", "SO2_a", "CO2_a", "NH3_a", "HNO3_a"]}}')

    data = netcdf.netcdf_file(outfile,"r")

    #plot
    plot_henry(data, "open", output_folder = "../outputs")

    #cleanup
    data.close()
    subprocess.call(["rm", "test_chem_dsl.nc"])

if __name__ == '__main__':
    main()
