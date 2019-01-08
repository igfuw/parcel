# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "../../")

from scipy.io import netcdf
import numpy as np
#import pytest
import math
import subprocess
import ast
import copy

from parcel import parcel
from libcloudphxx import common as cm
import functions as fn
from chem_conditions import parcel_dict


def plot_spectra(data, output_folder = '', output_title = ''):
    import Gnuplot

    # from ncdf file attributes read out_bin parameters as a dictionary ...
    out_bin = ast.literal_eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
#    assert out_bin["radii"]["lnli"] == 'lin', "this plot should be used with logarithmic spacing of bins"

    # left bin edges
    rd = data.variables["radii_r_dry"][:]
    drd = data.variables["radii_dr_dry"][:]
    rw = data.variables["cloud_r_wet"][:]
    drw = data.variables["cloud_dr_wet"][:]
#    drd = rd.copy()
#    for i in range(rd.shape[0]-1):
#      drd[i] = rd[i+1] - rd[i]
#    drd[rd.shape[0]-1] = drd[rd.shape[0] - 2]


    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    #d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10)

    g = Gnuplot.Gnuplot()# persist=1)
    g('set terminal pdf')
    g('set samples 500')

    ymin = 1e-7
    ymax = 100
    xmin = 0
    xmax = 45


    g('set output "' + output_folder + output_title + '.pdf"')
    g('set logscale xy')
#    g('invsqrt2pi = 0.398942280401433')
    g('lognormalx(x, sig_g, mean_g) = 1 / (x * log(sig_g) * sqrt(2 * pi)) * exp (-(log(x/mean_g))**2 / (2 * log(sig_g)**2))')
#    g('lognormallnx(x, sig, mean) = 1 / (log(sig) * sqrt(2 * pi)) * exp (-(x - log(mean))**2 / (2 * log(sig)**2))')
    g('set xrange[5e-3:40]')
    g('set yrange[1:3e4]')
    g('set multiplot')
#    g('plot \'ccn_clean_081018a3_with_widths\' u ($2*1e6):(($3/1e6)/($5*1e6)) w l lc 3')
   # g('plot  125 *  lognormalx(x, 1.20, 0.011) +  65 *  lognormalx(x, 1.7, 0.06) w l lc 2')
   # g('plot \'res_drop_cond\' u 3:($7/($4-$2)) w l lt 2 lc 6')
   # g('plot \'res_drop_cond_with_wet_radii\' u 12:($7/($13)) w l lt 3 lc 7')
    g('plot 125 *  lognormalx(x, 1.20, 0.011) +  65 *  lognormalx(x, 1.7, 0.06) w l lc 2, \'res_drop_cond\' u 3:($7/($4-$2)) w l lt 2 lc 6, \'res_drop_cond_with_wet_radii\' u 12:($7/($13)) w l lt 3 lc 7')# lt 3 lc 3 t \'dry radius\'')
#    g('plot \'res_drop_cond\' u 3:($5/200) w l lc 6')
#    g('plot 1')

#    for t in range(data.variables['t'].shape[0]):
    #g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
    #g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
    #g('set grid')
    #g('set nokey')
  
    #nd = data.variables['radii_m0'][t,:] * data.variables["rhod"][0] / 1.#d_log_rd
    nd = data.variables['radii_m0'][0,:] * data.variables["rhod"][0] /  drd
    nw = data.variables['cloud_m0'][0,:] * data.variables["rhod"][0] /  drw
    nw_end = data.variables['cloud_m0'][1,:] * data.variables["rhod"][0] /  drw
#    print nd
  
    plot_rd = Gnuplot.PlotItems.Data((rd+0.5*drd) * 1e6, nd * 1e-12, with_="l lt 3 lc 3", title="dry radius")
    plot_rw = Gnuplot.PlotItems.Data((rw+0.5*drw) * 1e6, nw * 1e-12, with_="line lw 1 lc 4", title="wet radius")
    plot_rw_end = Gnuplot.PlotItems.Data((rw+0.5*drw) * 1e6, nw_end * 1e-12, with_="line lw 1 lc 5", title="wet radius")
#    print plot_rd

#      print 'aaa'
#    g('show term')
    g.plot(plot_rd, plot_rw)
#    g.plot(plot_rw)
#    g.plot(plot_rw_end)
#      print 'aaaa'
    g('unset multiplot')


def main():

    data = netcdf.netcdf_file('Jon_GCCN.nc',   "r")

    # do the plotting
    plot_spectra(data, output_folder = "../outputs", output_title = "/Jon_GCCN_spectra")

    # cleanup
#    subprocess.call(["rm", p_dict['outfile']])

if __name__ == '__main__':
    main()
