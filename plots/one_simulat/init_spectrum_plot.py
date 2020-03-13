# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../../")
sys.path.insert(0, "../")
sys.path.insert(0, "./")

sys.path.insert(0, "../../unit_test")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest
import copy
import ast

from parcel import parcel
from chem_conditions import parcel_dict
import functions as fn

def plot_init_spectrum(data, outfolder):
    """
    Plot the initial dry diameter distribution and compare it with the analitycal solution

    """
    import Gnuplot

    # size distribution parameters from Kreidenweis 2003
    n_tot   = 125e6#, 15e6]#566e6
    mean_r  = 0.011e-6#, 0.14e-6]#0.04e-6
    gstdev  = 1.2#, 1.75]#2
    n_tot2   = 65e6#566e6
    mean_r2  = 0.06e-6#0.04e-6
    gstdev2  = 1.7#2

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = ast.literal_eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
    assert out_bin["drad"]["lnli"] == 'log', "this test should be run with logarithmic spacing of bins"

    # parcel initial condition
    rd = data.variables["drad_r_dry"][:] # left bin edges

    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10)

    # initial size distribution from the model
    model = data.variables['drad_m0'][0,:] * data.variables["rhod"][0] / d_log_rd

    # variables for plotting theoretical solution
    radii = np.logspace(-3, 1, 100) * 1e-6
    theor = np.empty(radii.shape)
    theor2 = np.empty(radii.shape)
    for it in range(radii.shape[0]):
        theor[it] = fn.log10_size_of_lnr(n_tot, mean_r, math.log(radii[it], 10), gstdev)
        theor2[it] = fn.log10_size_of_lnr(n_tot2, mean_r2, math.log(radii[it], 10), gstdev2)
    g = Gnuplot.Gnuplot()
    g('set term svg dynamic enhanced')
    g('reset')
    g('set output "' + outfolder + '/init_spectrum.svg" ')
    g('set logscale x')
    g('set xlabel "particle dry diameter [μm]" ')
    g('set ylabel "dN/dlog_{10}(D) [cm^{-3} log_{10}(size interval)]"')
    g('set grid')
    g('set xrange [0.001:10]')
    g('set yrange [0:800]')

    theory_r = Gnuplot.PlotItems.Data(radii * 2 * 1e6,  theor * 1e-6, with_="lines", title="theory")
    theory_r2 = Gnuplot.PlotItems.Data(radii * 2 * 1e6,  theor2 * 1e-6, with_="lines", title="theory2")
    plot     = Gnuplot.PlotItems.Data(rd    * 2 * 1e6,  model * 1e-6, with_="steps", title="model" )

    g.plot(theory_r, theory_r2, plot)

def main():
    # copy options from chem_conditions ...
    opts_dict = copy.deepcopy(parcel_dict)

    # ... and modify them for the current test
    opts_dict['outfile']  = "test_init_spectrum.nc"
    opts_dict['chem_dsl'] = True

    opts_dict['z_max']    = .05
    opts_dict['dt']       = .1
    opts_dict['w']        = .5
    opts_dict['outfreq']  = 1
    opts_dict['sd_conc']  = 1024 * 44
    opts_dict['outfreq']  = 1

    opts_dict['out_bin'] = '{"drad": {"rght": 2.5e-06, "left": 5e-09, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [0,1,3]}}'

                            #{"rght": 2.5e-05, "moms": [0,1,3], "drwt": "wet", "nbin": 26, "lnli": "log", "left": 5e-07}}

    # run parcel
    parcel(**opts_dict)

    # simulation results
    data = netcdf.netcdf_file(opts_dict['outfile'],   "r")

    # doing plotting 
    plot_init_spectrum(data, outfolder="../outputs/")

    # cleanup
    subprocess.call(["rm", opts_dict['outfile']])

if __name__ == '__main__':
    main()
