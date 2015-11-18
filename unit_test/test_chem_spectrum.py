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
    opts_dict['sd_conc']  = 1024

    opts_dict['out_bin'] = '{"drad": {"rght": 1e-06, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [0]}}'

    # run parcel
    #parcel(**opts_dict)

    # simulation results
    data = netcdf.netcdf_file(opts_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", opts_dict['outfile']])

    #request.addfinalizer(removing_files)
    return data

def test_init_spectrum(data, eps = 1):
    """
    Checking if the initial aerosol size distribution is close to the analytic solution 

    """
    import Gnuplot

    ymax = 1e4
    ymin = 1e-1

    # log size distribution
    def log_size_of_r(n_tot, mean_r, radii, gstdev):
        return n_tot / math.sqrt(2.* math.pi) / math.log(gstdev, 10) / radii\
               * math.exp(-1. * math.pow(math.log(radii/mean_r, 10) / math.sqrt(2.) / math.log(gstdev, 10), 2))

    def log_size_of_lnr(n_tot, mean_r, lnr, gstdev):
        return n_tot / math.sqrt(2 * math.pi) / math.log(gstdev, 10)\
               * math.exp(-1. * math.pow((lnr - math.log(mean_r, 10)), 2) / 2. / math.pow(math.log(gstdev,10),2))

    # dry diameter bins
    rd      = data.variables["drad_r_dry"][:] * 1e6
    drd     = np.empty(rd.shape)
    drd[0]  = rd[0] - 0
    drd[1:] = (rd[1:] - rd[0:-1])

    # initial dry radii
    nd = data.variables['drad_m0'][0,:] / drd / 1e6

    # theoretical dry radii
    td  = np.empty(rd.shape)
    td2 = np.empty(rd.shape)

    for i in range(td.shape[0]):
      td[i]  = log_size_of_r(  n_tot, mean_r, rd[i] / 1e6, gstdev) / drd[i] / 1e6
      #td2[i] = log_size_of_lnr(n_tot, mean_r, rd[i] / 1e6, gstdev) / drd[i] / 1e6

    g = Gnuplot.Gnuplot()
    g('set term svg dynamic enhanced')
    g('reset')
    g('set output "plots/outputs/init_spectrum.svg" ')
    g('set logscale xy')
    g('set ylabel "[mg^{-1} μm^{-1}]"')
    #g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
    g('set grid')
    g('set nokey')

    g('set xlabel "particle dry diameter [μm]" ')

    plot_rd  = Gnuplot.PlotItems.Data(rd, nd, with_="fsteps", title="dry diameter")
    plot_td  = Gnuplot.PlotItems.Data(rd, td, with_="fsteps", title="theory")
    plot_td2 = Gnuplot.PlotItems.Data(rd, td2, with_="fsteps", title="theory 2")

    g.plot(plot_rd, plot_td, plot_td2)

