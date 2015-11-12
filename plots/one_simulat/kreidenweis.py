# This Python file uses the following encoding: utf-8
import sys
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
from functions import *

def plot_fig1(data, output_folder = '', output_title = ''):

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    spn_idx = 0
    #spn_idx = int(math.ceil(float(f_out_chem['open'].chem_spn)/float(f_out_chem['open'].outfreq)))

    # plot settings
    plt.figure(1, figsize=(10,12))
    plots = []
    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)
    for ax in plots:
      ax.set_ylabel('t [s]')

    # read in y-axis (time)
    t   = data.variables["t"][spn_idx:]

    # calculate lwc
    plots[0].set_xlabel('lwc g/kg dry air')
    plots[0].grid()
    plots[0].set_xlim([0., 2.5])
    plots[0].set_xticks([0., 0.5, 1, 1.5, 2, 2.5])
    plots[0].plot(np.sum(data.variables["radii_m3"][spn_idx:], axis=1) * 4. / 3 * math.pi * 998.2 * 1e3, t, "b.-")

    # calculate SO2 gas volume concentration
    p    = data.variables["p"][spn_idx:]
    T    = data.variables["T"][spn_idx:]
    rhod = data.variables["rhod"][spn_idx:]

    plots[1].set_xlabel('SO2 conc (ppb) - TODO aq')
    plots[1].grid()
    plots[1].set_xlabel('gas vol.conc SO2 [ppb]')
    plots[1].set_xticks([0., 0.05, 0.1, 0.15, 0.2])
    plots[1].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2'])
    plots[1].set_xlim([0., 0.2])
    tmp1 = mix_ratio_to_mole_frac(data.variables["SO2_g"][spn_idx:], p, cm.M_SO2,     T, rhod)
    tmp2 = mix_ratio_to_mole_frac(data.variables["SO2_a"][spn_idx:], p, cm.M_SO2_H2O, T, rhod)
    plots[1].plot(tmp1 * 1e9, t, "r.-")
    plots[1].plot(tmp2 * 1e9, t, "b.-")
    plots[1].plot((tmp1 + tmp2) * 1e9, t, "g.-")

    print data.variables["SO2_a"][spn_idx:]

    # calculate average pH
    # (weighted with volume of cloud droplets)
    plots[2].set_xlabel('average pH')
    plots[2].grid()

    r3     = data.variables["radii_m3"][spn_idx:]
    n_H    = data.variables["chem_H"][spn_idx:] / cm.M_H
    nom    = np.zeros(t.shape[0])
    den    = np.zeros(t.shape[0])
    for time in range(t.shape[0]): 
        for it, val in enumerate(r3[time,:]):
            if val > 0:
                 nom[time] += (n_H[time, it] / (4./3 * math.pi * val * 1e3)) * val
        den[time] = np.sum(r3[time,:])                                 # to liters

    pH  = -1 * np.log10(nom / den)
    # plots[2].set_xlim([3.6, 4.8])
    # plots[2].set_xticks([3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8])
    plots[2].plot(pH, t, "b.-")

    plt.savefig(output_folder + output_title + ".pdf")

def plot_fig2(data, output_folder = '', output_title = ''):
    import Gnuplot

    g = Gnuplot.Gnuplot()# persist=1)
    g('set term svg dynamic enhanced')

    ymax = 1e4
    ymin = 1e-1

    # diameter
    rw = data.variables["radii_r_wet"][:] * 1e6 * 2  
    rd = data.variables["radiidry_r_dry"][:] * 1e6 * 2

    #TODO - add test if it is == to dr in netcdf
    drw = np.empty(rw.shape)
    drw[0] = rw[0] - 0
    drw[1:] = (rw[1:] - rw[0:-1]) * 1e6

    drd = np.empty(rd.shape)
    drd[0] = rd[0] - 0
    drd[1:] = (rd[1:] - rd[0:-1]) * 1e6

    for t in range(data.variables['t'].shape[0]):

        g('reset')
        g('set output "' + output_folder + '/Kreidenweis_plot_spec_' + str("%03d" % t) + '.svg"')
        g('set logscale xy')
        g('set ylabel "[mg^{-1} μm^{-1}]"')
        g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
        g('set grid')
        g('set nokey')

        # FSSP range
        g('set arrow from .5,' + str(ymin) + 'to .5,' + str(ymax) + 'nohead')
        g('set arrow from 25,' + str(ymin) + 'to 25,' + str(ymax) + 'nohead')

        g('set xlabel "particle diameter [μm]" ')

        nw = data.variables['radii_m0'][t,:] / drw
        nd = data.variables['radiidry_m0'][t,:] / drd

        plot_rw = Gnuplot.PlotItems.Data(rw, nw, with_="fsteps", title="wet radius")
        plot_rd = Gnuplot.PlotItems.Data(rd, nd, with_="fsteps", title="dry radius")

        g.plot(plot_rw, plot_rd)
 
def main():
    # initial condition
    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = rh_to_rv(RH_init, T_init, p_init)

    # calculate rhod for initial gas mixing ratio
    rhod_init   = rhod_calc(T_init, p_init, r_init)
    # initial condition for trace geses
    SO2_g_init  = mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
    O3_g_init   = mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
    H2O2_g_init = mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
    CO2_g_init  = mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
    NH3_g_init  = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
    HNO3_g_init = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

    # aerosol size distribution
    mean_r = .08e-6 / 2
    gstdev = 2.
    n_tot  = 566.e6

    # process toggling
    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    # output
    z_max       = 2400.
    dt          = .1
    w           = 1.
    outfreq     = int(z_max / dt / 100)
    sd_conc     = 1024.
    outfile     = "Kreidenweis.nc"

    # define output for moments and chemistry
    out_bin = '{"radii":    {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500, "moms": [0, 3]},\
                "radiidry": {"rght": 1e-4, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 500, "moms": [0, 3]},\
                "chem" :    {"rght": 1e-4, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 500,\
                                   "moms": ["H", "OH", "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI"]}}'

    # run parcel, run!
    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
           chem_sys = 'closed', outfile = outfile,\
           sd_conc = sd_conc,\
           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
           out_bin = out_bin)

    # TODO - why do I have to repeat this import here?
    from scipy.io import netcdf
    data = netcdf.netcdf_file(outfile,   "r")

    plot_fig1(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig1")
    plot_fig2(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig2")

    #cleanup
    #subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()
