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

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    # plot settings
    plt.figure(1)
    plt.rcParams.update({'font.size': 8})
    plots = []
    for i in range(3):
      plots.append(plt.subplot(1,3,i+1))
                             #(rows, columns, number)
    for ax in plots:
      ax.set_ylabel('t [s]')
      ax.set_ylim([0, 2400])
      ax.set_yticks([0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400])

    #plt.tight_layout()
    spn_idx = 0#4

    # read in y-axis (time)
    t   = data.variables["t"][spn_idx:] - data.variables["t"][spn_idx]

    # calculate lwc
    plots[0].set_xlabel('lwc g/kg dry air')
    plots[0].grid()
    #plots[0].set_xlim([0., 2.5])
    #plots[0].set_xticks([0., 0.5, 1, 1.5, 2, 2.5])
    plots[0].plot(np.sum(data.variables["radii_m3"][spn_idx:], axis=1) * 4. / 3 * math.pi * 998.2 * 1e3, t, "b.-")

    # calculate SO2 gas volume concentration
    p    = data.variables["p"][spn_idx:]
    T    = data.variables["T"][spn_idx:]
    rhod = data.variables["rhod"][spn_idx:]

    plots[1].grid()
    plots[1].set_xlabel('gas vol.conc SO2 [ppb]')
    plots[1].set_xticks([0., 0.05, 0.1, 0.15, 0.2])
    plots[1].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2'])
    plots[1].set_xlim([0., 0.2])
    tmp1 = fn.mix_ratio_to_mole_frac(data.variables["SO2_g"][spn_idx:], p, cm.M_SO2, T, rhod)
    tmp2 = fn.mix_ratio_to_mole_frac(\
      np.squeeze(data.variables["plt_ch_SO2_a"][spn_idx:]), p, cm.M_SO2_H2O, T, rhod)
    #tmp2 = fn.mix_ratio_to_mole_frac(data.variables["SO2_a"][spn_idx:], p, cm.M_SO2_H2O, T, rhod)
    plots[1].plot((tmp1 + tmp2) * 1e9, t, "g.-")
    #plots[1].plot((tmp2) * 1e9, t, "r.-")
    #plots[1].plot((tmp1) * 1e9, t, "b.-")

    # calculate average pH
    # (weighted with volume of cloud droplets)
    plots[2].set_xlabel('average pH')
    plots[2].grid()
    #plots[2].set_xlim([3.8, 5])
    #plots[2].set_xticks([3.8, 4, 4.2, 4.4, 4.6, 4.8, 5])
 
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
    plots[2].plot(pH, t, "b.-")

    plt.savefig(output_folder + output_title + ".pdf")

def plot_fig2(data, output_folder = '', output_title = ''):
    import Gnuplot

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = ast.literal_eval(getattr(data, "out_bin"))
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
    ymax = 1500
    xmin = 0.001
    xmax = 10

    for t in range(data.variables['t'].shape[0]):
        if t % 10 == 0:
            g('reset')
            g('set output "' + output_folder + '/Kreidenweis_plot_spec_' + str("%03d" % t) + '.svg"')
            g('set logscale x')
            g('set xlabel "particle diameter [μm]" ')
            g('set ylabel "dN/dlog_{10}(D) [cm^{-3} log_{10}(size interval)]"')
            g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
            g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
            g('set grid')
            g('set nokey')
    
            nd = data.variables['specd_m0'][t,:] * data.variables["rhod"][0] / d_log_rd
    
            plot_rd = Gnuplot.PlotItems.Data(rd * 2 * 1e6, nd * 1e-6, with_="steps", title="dry radius")
    
            g.plot(plot_rd)

def plot_fig3(data, output_folder = '', output_title = ''):
    import Gnuplot

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = ast.literal_eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
    assert out_bin["specd"]["lnli"] == 'log', "this plot should be used with logarithmic spacing of bins"

    # left bin edges
    rd = data.variables["specd_r_dry"][:]

    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10)

    g = Gnuplot.Gnuplot()# persist=1)
    g('set term svg dynamic enhanced')

    ymin = .01
    ymax = 10
    xmin = 0.01
    xmax = 1

    g('reset')
    g('set output "' + output_folder + output_title + '.svg"')
    g('set logscale xy')
    g('set xlabel "particle diameter [μm]" ')
    g('set ylabel "dS(VI)/dlog_{10}(D) [μg /m^3 log_{10}(size interval)]"')
    g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
    g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
    g('set grid')
    g('set nokey')

    s6_ini = data.variables['chemd_S_VI'][0,:]  * data.variables["rhod"][0]  / d_log_rd
    s6_end = data.variables['chemd_S_VI'][-1,:] * data.variables["rhod"][-1] / d_log_rd

    mass1 = 1.35837239e-10 * 1e9
    mass2 = 2.00500078e-09 * 1e9
    edge1 = .08
    edge2 = .125

    #g('set arrow from ' + str(edge1) + ',' + str(mass1) + 'to ' + str(edge1) + ',' + str(ymin) + 'nohead lw 2')
    #g('set arrow from ' + str(edge2) + ',' + str(mass2) + 'to ' + str(edge2) + ',' + str(ymin) + 'nohead lw 2')
 
    plot_ini = Gnuplot.PlotItems.Data(rd * 2 * 1e6, s6_ini * 1e9, with_="steps lw 3 lt 1", title="ini")
    plot_end = Gnuplot.PlotItems.Data(rd * 2 * 1e6, s6_end * 1e9, with_="steps lw 3 lt 2", title="end")

    g.plot(plot_ini, plot_end)

def plot_fig4(data, output_folder = '', output_title = ''):

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    # plot settings
    plt.figure(1, figsize=(18,14))
    plt.rcParams.update({'font.size': 12})
    plots = []
    for i in range(6):
      plots.append(plt.subplot(2,3,i+1))
                             #(rows, columns, number)
    for ax in plots:
      ax.set_ylabel('t [s]')
      ax.set_ylim([0, 2400])
      ax.set_yticks([0, 400, 800, 1200, 1600, 2000, 2400])

    #plt.tight_layout()
    spn_idx = 4

    # read in y-axis (time)
    t   = data.variables["t"][spn_idx:] - data.variables["t"][spn_idx]

    # calculate SO2 gas volume concentration
    p    = data.variables["p"][spn_idx:]
    T    = data.variables["T"][spn_idx:]
    rhod = data.variables["rhod"][spn_idx:]

 
    chem_dict = {"SO2" : {"it": 0, "unit": "[ppb]", "scale": 1e9, "lim": [0.06, 0.2],   "Mg": cm.M_SO2,  "Ma": cm.M_SO2_H2O},\
                 "O3"  : {"it": 1, "unit": "[ppb]", "scale": 1e9, "lim": [50., 50.1],   "Mg": cm.M_O3,   "Ma": cm.M_O3},\
                 "H2O2": {"it": 2, "unit": "[ppb]", "scale": 1e9, "lim": [0.45, 0.5],   "Mg": cm.M_H2O2, "Ma": cm.M_H2O2},\
                 "CO2" : {"it": 3, "unit": "[ppm]", "scale": 1e6, "lim": [360., 361.2], "Mg": cm.M_CO2,  "Ma": cm.M_CO2_H2O},\
                 "NH3" : {"it": 4, "unit": "[ppb]", "scale": 1e9, "lim": [0., 0.1],     "Mg": cm.M_NH3,  "Ma": cm.M_NH3_H2O},\
                 "HNO3": {"it": 5, "unit": "[ppb]", "scale": 1e9, "lim": [0., 0.1],     "Mg": cm.M_HNO3, "Ma": cm.M_HNO3}}

    for key, val in chem_dict.iteritems():

        plots[val["it"]].grid()
        plots[val["it"]].set_xlabel('gas + aq vol.conc' + key + val["unit"])
        #plots[1].set_xticks([0., 0.05, 0.1, 0.15, 0.2])
        #plots[1].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2'])
        plots[val["it"]].set_xlim(val["lim"])
        tmp1 = fn.mix_ratio_to_mole_frac(data.variables[key + "_g"][spn_idx:], p, val["Mg"], T, rhod)
        tmp2 = fn.mix_ratio_to_mole_frac(data.variables[key + "_a"][spn_idx:], p, val["Ma"], T, rhod)
        plots[val["it"]].plot((tmp1 + tmp2) * val["scale"], t, "g.-")

    plt.savefig(output_folder + output_title + ".pdf")

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
                  "radii" : {"rght": 1e-4, "left": 1e-10, "drwt": "wet", "lnli": "log", "nbin": 100, "moms": [3]},\
                   "specd": {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 100, "moms": [0, 1, 3]},\
                  "plt_rw": {"rght": 1,    "left": 0,     "drwt": "wet", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_rd": {"rght": 1,    "left": 0,     "drwt": "dry", "lnli": "lin", "nbin": 1,   "moms": [0, 1, 3]},\
                  "plt_ch": {"rght": 1,    "left": 0,     "drwt": "dry", "lnli": "lin", "nbin": 1,\
                             "moms": ["O3_a", "H2O2_a", "H", "OH", "SO2_a", "S_VI", "CO2_a", "NH3_a", "HNO3_a"]}}'

    # run parcel
    #parcel(**p_dict)

    # simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # do the plotting
    plot_fig1(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig1")
    plot_fig2(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig2")
    plot_fig3(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig3")
    plot_fig4(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig4")

    # cleanup
    subprocess.call(["rm", p_dict['outfile']])

if __name__ == '__main__':
    main()
