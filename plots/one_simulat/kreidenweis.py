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
import ast
import copy

from parcel import parcel
from libcloudphxx import common as cm
import functions as fn
from chem_conditions import parcel_dict

def plot_fig1(data, output_folder = '', output_title = ''):

    # read in data
    spn_idx = 2 * 5
    t    = data.variables["t"][:] #- data.variables["t"][spn_idx]
    p    = data.variables["p"][:]
    T    = data.variables["T"][:]
    rhod = data.variables["rhod"][:]
    r_v  = data.variables["r_v"][:]
    rho  = fn.rho_calc(T, p, r_v)

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    # plot settings
    fig = plt.figure(figsize=(28,13))
    plt.rcParams.update({'font.size': 30})
    y_labels  = [0,   200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400] #time above CB
    y_labels2 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6,   0.7, 0.8,   0.9,    1, 1.1,  1.2,  1.3 ] #height above CB

    # calculate lwc
    ax = fig.add_subplot(131)
    ax.set_xlim([0., 2.5])
    ax.set_ylim([200, 2600])
    ax.set_yticks([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600])
    ax.set_yticklabels(y_labels)

    ax.tick_params(axis='x', pad=15)
    ax.tick_params(axis='y', pad=15)

    ax.set_xlabel('LWC [g/kg]')
    ax.set_ylabel('time above cloud base [s]')

    plt.grid()
    plt.plot(data.variables["acti_m3"][:] * 4./3 * math.pi * 999.5 * 1e3 * rhod / rho, t,"b", lw=4.)

    # calculate SO2 gas volume concentration
    tmp1 = fn.mix_ratio_to_mole_frac(data.variables["SO2_g"][:], p, cm.M_SO2, T, rhod)
    tmp2 = fn.mix_ratio_to_mole_frac(\
      np.squeeze(data.variables["plt_ch_SO2_a"][:]), p, cm.M_SO2_H2O, T, rhod)

    ax = fig.add_subplot(132)
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax.set_ylim([200, 2600])
    ax.set_yticks([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600])
 
    ax.set_xlim([0., 0.2])
    ax.set_xticks([0., 0.05, 0.1, 0.15, 0.2])
    ax.set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2'])

    ax.set_xlabel('$\mathrm{SO_2}$ conc. [ppb]') #gas volume concentration

    plt.grid()
    plt.plot((tmp1 + tmp2) * 1e9 * rhod / rho, t, "b", lw=4.)
    print("total % of SO2 that was converted to H2SO4 = ", (1 - (tmp1[-1] + tmp2[-1]) / (tmp1[0] + tmp2[0])) * 100)

    # calculate average pH
    # (weighted with volume of cloud droplets)
    r3     = data.variables["radii_m3"][:]
    n_H    = data.variables["chem_H"][:] / cm.M_H
    nom    = np.zeros(t.shape[0])
    den    = np.zeros(t.shape[0])
    for time in range(t.shape[0]): 
        for it, val in enumerate(r3[time,:]):
            if val > 0:
                 nom[time] += (n_H[time, it] / (4./3 * math.pi * val * 1e3)) * val
        den[time] = np.sum(r3[time,:])                                 # to liters

    pH  = -1 * np.log10(nom / den)

    ax = fig.add_subplot(133)

    ax.set_ylim([200, 2600])
    ax.set_yticks([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600])
    ax.set_yticklabels(y_labels2)
    ax.yaxis.tick_right()

    ax.set_xlabel('average pH')
    ax.set_ylabel('height above cloud base [km]')
    ax.yaxis.set_label_position("right")

    ax.set_xlim([3.6, 5])
    ax.set_xticks([3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5])

    plt.grid()
    plt.plot(pH, t, "b", lw=4.)

    plt.savefig(output_folder + output_title + ".pdf")

def plot_fig2(data, output_folder = '', output_title = ''):
    import Gnuplot

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
    assert out_bin["specd"]["lnli"] == 'log', "this plot should be used with logarithmic spacing of bins"

    # left bin edges
    rd = data.variables["specd_r_dry"][:]

    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10)

    g = Gnuplot.Gnuplot()# persist=1)
    g('set term svg dynamic enhanced')

    ymin = 0
    ymax = 1200
    xmin = 0.01
    xmax = 1

    for t in range(data.variables['t'].shape[0]):
        if t % 10 == 0:
            g('reset')
            g('set output "' + output_folder + output_title +'_plot_spec_' + str("%03d" % t) + '.svg"')
            g('set logscale x')
            g('set xlabel "particle diameter [μm]" ')
            g('set ylabel "dN/dlog_{10}(D) [cm^{-3} log_{10}(size interval)]"')
            g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
            g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
            g('set grid')
            g('set nokey')
    
            nd = data.variables['specd_m0'][t,:] * data.variables["rhod"][0] / d_log_rd
    
            plot_rd = Gnuplot.PlotItems.Data(rd * 2 * 1e6, nd * 1e-6, with_="steps lw 2", title="dry radius")
    
            g.plot(plot_rd)

def plot_fig3(data, output_folder = '', output_title = ''):
    import Gnuplot

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
    assert out_bin["specd"]["lnli"] == 'log', "this plot should be used with logarithmic spacing of bins"

    # left bin edges
    rd = data.variables["specd_r_dry"][:]

    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10)

    g = Gnuplot.Gnuplot()# persist=1)

    s6_ini = data.variables['chemd_S_VI'][0,:]  * data.variables["rhod"][0]  / d_log_rd
    s6_end = data.variables['chemd_S_VI'][-1,:] * data.variables["rhod"][-1] / d_log_rd

    plot_ini = Gnuplot.PlotItems.Data(rd * 2 * 1e6, s6_ini * 1e9, with_="steps lw 3 lt 1", title="ini")
    plot_end = Gnuplot.PlotItems.Data(rd * 2 * 1e6, s6_end * 1e9, with_="steps lw 3 lt 2", title="end")

    xmin = 0.01
    xmax = 1
    ymin = .01
    ymax = 15

    g('reset')
    g('set term svg dynamic enhanced font "Verdana, 14"')
    g('set xlabel "particle diameter [μm]" ')
    g('set ylabel "dS(VI)/dlog_{10}(D) [μg /m^3 per log_{10}(size interval)]"')
    g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
    g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
    g('set grid')
    g('set nokey')

    g('set output "' + output_folder + output_title + '.svg"')
    g('set logscale xy')

    g.plot(plot_ini, plot_end)

    g('reset')
    g('set term svg dynamic enhanced font "Verdana, 14"')
    g('set xlabel "particle diameter [μm]" ')
    g('set ylabel "dS(VI)/dlog_{10}(D) [μg /m^3 per log_{10}(size interval)]"')
    g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
    g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
    g('set ytics 0,5,15')
    g('set grid')
    g('set nokey')

    g('set output "' + output_folder + output_title + '_no_ylog.svg"')
    g('set logscale x')

    g.plot(plot_ini, plot_end)


def plot_pH_size_dist(data, output_folder = '', output_title = ''):
    import Gnuplot

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
    assert out_bin["chem"]["lnli"] == 'log', "this plot should be used with logarithmic spacing of bins"

    # left bin edges
    rw = data.variables["chem_r_wet"][:]

    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    d_log_rw = math.log(rw[2], 10) - math.log(rw[1], 10)


    r3     = data.variables["radii_m3"][-1]
    n_H    = data.variables["chem_H"][-1] / cm.M_H
    conc_H = np.ones(r3.shape)
    for it, val in enumerate(r3[:]):
        if val > 0:
             conc_H[it] = n_H[it] / (4./3 * math.pi * val * 1e3)

    pH  = -1 * np.log10(conc_H)

    g = Gnuplot.Gnuplot()# persist=1)
    g('set term svg dynamic enhanced font "Verdana, 14"')

    ymin = 0
    ymax = 6
    xmin = 10
    xmax = 100

    g('reset')
    g('set output "' + output_folder + output_title + '.svg"')
    g('set logscale x')
    g('set xlabel "particle diameter [μm]" ')
    g('set ylabel "pH"')
    g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
    g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
    g('set grid')
    g('set nokey')

    plot_end = Gnuplot.PlotItems.Data(rw * 2 * 1e6, pH , with_="steps lw 3 lt 2")

    g.plot(plot_end)

def main():

    # copy options from chem_conditions ...
    p_dict = copy.deepcopy(parcel_dict)

    # ... and modify them for the current test
    p_dict['outfile']  = "Kreidenweis.nc"

    p_dict['chem_dsl'] = True
    p_dict['chem_dsc'] = True
    p_dict['chem_rct'] = True

    p_dict['sd_conc']  = 1025
    p_dict['outfreq']  = int(p_dict['z_max'] / p_dict['dt'] / 100) * 4

    p_dict['out_bin']  = '{\
                  "chem"  : {"rght": 1e-4, "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 100, "moms": ["H"]},\
                  "chemd" : {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 100, "moms": ["S_VI"]},\
                  "radii" : {"rght": 1e-4, "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 100, "moms": [0, 1, 3]},\
                   "specd": {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 100, "moms": [0, 1, 3]},\
                  "plt_rw": {"rght": 1,    "left": 0,     "drwt": "wet", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_rd": {"rght": 1,    "left": 0,     "drwt": "dry", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_ch": {"rght": 1,    "left": 0,     "drwt": "dry", "lnli": "lin", "nbin": 1,\
                             "moms": ["O3_a", "H2O2_a", "H", "OH", "SO2_a", "S_VI", "CO2_a", "NH3_a", "HNO3_a"]}}'

    # run parcel
    parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # do the plotting
    plot_fig1(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig1")
    plot_fig2(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig2")
    plot_fig3(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig3")

    # cleanup
    subprocess.call(["rm", p_dict['outfile']])

if __name__ == '__main__':
    main()
