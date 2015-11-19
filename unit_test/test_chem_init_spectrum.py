# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest
import copy

from parcel import parcel
from libcloudphxx import common as cm
from chemical_plot import plot_chem
from chem_conditions import *

@pytest.fixture(scope="module")
def data(request):

    # copy options from chem_conditions ...
    opts_dict = copy.deepcopy(parcel_dict)

    # ... and modify them for the current test
    opts_dict['outfile']  = "test_chem_spectrum.nc"
    opts_dict['chem_dsl'] = True

    opts_dict['z_max']    = .5
    opts_dict['dt']       = .1
    opts_dict['w']        = .5
    opts_dict['outfreq']  = 1
    opts_dict['sd_conc']  = 1024 * 50

    opts_dict['out_bin'] = '{"drad": {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 50, "moms": [0]},\
                           "wradii": {"rght": 1e-4, "left": 1e-10, "drwt": "wet", "lnli": "lin", "nbin": 50, "moms": [0, 3]},\
                           "dradii": {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "lin", "nbin": 50, "moms": [0, 3]}}'

    # run parcel
    parcel(**opts_dict)

    # simulation results
    data = netcdf.netcdf_file(opts_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", opts_dict['outfile']])

    #request.addfinalizer(removing_files)
    return data

def test_todo(data, eps = 1e-20):
    """
    Check if the initial mass is the same as in the article
    TODO - convergence with the initial sd_conc

    """
    chem_rho = getattr(data, "chem_rho")
    mom3     = data.variables["dradii_m3"][:,:]
    rhod     = data.variables["rhod"][:]
    rv       = data.variables["r_v"][:]

    ini = mom3[0,:].sum()  * 4./3 * math.pi * chem_rho
    end = mom3[-1,:].sum() * 4./3 * math.pi * chem_rho

    rhod_parc_init = rhod[0]

    print " "

    print "sd_conc = 1024 * 50,   nbin = 50,  mix = 2.04312254645e-09  [kg/kg dry air]"
    print "sd_conc = 1024 * 100,  nbin = 50,  mix = 2.02454798163e-09  [kg/kg dry air]"
    print "sd_conc = 1024 * 200,  nbin = 50,  mix = 2.00070482127e-09  [kg/kg dry air]"
    print "sd_conc = 1024 * 300,  nbin = 50,  mix = 1.98218151415e-09  [kg/kg dry air]"
    print "sd_conc = 1024 * 500,  nbin = 50,  mix = 1.97226250868e-09  [kg/kg dry air]"
    print "sd_conc = 1024 * 1000, nbin = 50,  mix = 1.92057052417e-09  [kg/kg dry air]"
    print "sd_conc = 1024 * 2000, nbin = 50,  mix = 1.87114810921e-09  [kg/kg dry air]"

    print " "
    print "initial mixing ratio: ", ini, " [kg/kg dry air]"
    print "final mixing ratio:   ", end, " [kg/kg dry air]"
    print " "
 
    a = 2.375 * 1e-9 / rhod_parc_init
    print "article init mixing ratio (/ rhod parc init):  ", a, " [kg / kg dry air]"

def test_plot_init_spectrum(data, eps = 1):
    """
    Checking if the initial aerosol size distribution is close to the analytic solution 

    """
    import Gnuplot

    # size distribution params
    n_tot   = 566e6
    mean_r  = 0.04e-6
    gstdev  = 2

    # size distribution (defined as a function of log_10(radius))
    def log10_size_of_lnr(n_tot, mean_r, lnr, gstdev):
        return n_tot / math.sqrt(2 * math.pi) / math.log(gstdev, 10)\
               * math.exp(-1. * math.pow((lnr - math.log(mean_r, 10)), 2) / 2. / math.pow(math.log(gstdev,10),2))

    # variables for plotting theoretical solution
    # TODO - use the same bins as in the model
    radii     = np.logspace(-3, 1, 100) * 1e-6
    radii_log = np.empty(radii.shape)
    theor_n     = np.empty(radii.shape) 
    for it in range(radii.shape[0]):
        radii_log[it] = math.log(radii[it], 10)
        theor_n[it]     = log10_size_of_lnr(n_tot, mean_r, radii_log[it], gstdev)

    # variables for plotting model initail condition
    # (for plotting purpose model solution needs to be divided by ln(d1) - ln(d2) )
    rd       = data.variables["drad_r_dry"][:]
    log_rd   = np.empty(rd.shape)
    d_log_rd = np.empty(rd.shape)

    for it in range(rd.shape[0]):
        log_rd[it] = math.log(rd[it], 10)

    d_log_rd[0]  = log_rd[0] - 1
    d_log_rd[1:] = (log_rd[1:] - log_rd[0:-1])

    # initial size distribution from the model
    model_n = data.variables['drad_m0'][0,:]

    # initial density
    rhod = data.variables["rhod"][0]

    g = Gnuplot.Gnuplot()
    g('set term svg dynamic enhanced')
    g('reset')
    g('set output "plots/outputs/init_spectrum.svg" ')
    g('set logscale x')
    g('set xlabel "particle dry diameter [μm]" ')
    g('set ylabel "dN/dlog_{10}(D) [cm^{-3} log_{10}(size interval)]"')
    g('set grid')
    g('set xrange [0.001:10]')
    g('set yrange [0:800]')

    theory_r = Gnuplot.PlotItems.Data(radii * 2 * 1e6,  theor_n   * 1e-6,                 with_="lines", title="theory")
    plot     = Gnuplot.PlotItems.Data(rd    * 2 * 1e6,  model_n * rhod * 1e-6 / d_log_rd, with_="steps", title="model" )

    g.plot(theory_r, plot)


