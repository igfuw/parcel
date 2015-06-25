# This Python file uses the following encoding: utf-8

import sys
import Gnuplot
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
import matplotlib.pyplot as plt
from scipy.io import netcdf
from parcel import parcel
import numpy as np
import pytest

def test_spectrum_plot():
    # running parcel model for open / closed chem system  ...
    parcel(dt = 1., sd_conc = 512, outfreq = 10,  outfile="test_spectrum.nc")

    # ... plotting the results ...
    f_out = netcdf.netcdf_file("test_spectrum.nc", "r")

    g = Gnuplot.Gnuplot()# persist=1)
    g('set term svg dynamic enhanced')

    ymax = 1e4
    ymin = 1e-1

    r = f_out.radii * 1e6
    dr = np.empty(r.shape) 
    dr[0] = r[0] - 0
    dr[1:] = (r[1:] - r[0:-1]) * 1e6

    for t in range(f_out.variables['t'].shape[0]):

        g('reset')
        g('set output "plot_spec_' + str("%03d" % t) + '.svg"')
        g('set logscale xy')
        g('set ylabel "[mg^{-1} μm^{-1}]"')
        g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
        g('set grid')
        g('set nokey')

        # FSSP range
        g('set arrow from .5,' + str(ymin) + 'to .5,' + str(ymax) + 'nohead')
        g('set arrow from 25,' + str(ymin) + 'to 25,' + str(ymax) + 'nohead')

        g('set xlabel "particle radius [μm]" ')

        n = f_out.variables['conc'][t,:] / dr
       
        g.plot(Gnuplot.Data(r, n, with_='fsteps')) #... _bins

test_spectrum_plot()
