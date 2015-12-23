# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../../")
sys.path.insert(0, "../")
sys.path.insert(0, "./")

from scipy.io import netcdf
import numpy as np
import pytest
import subprocess

from parcel import parcel

def plot_spectrum(data, outfolder):
    import Gnuplot

    g = Gnuplot.Gnuplot()# persist=1)
    g('set term svg dynamic enhanced')

    ymax = 1e4
    ymin = 1e-1

    rw = data.variables["wradii_r_wet"][:] * 1e6
    rd = data.variables["dradii_r_dry"][:] * 1e6

    #TODO - add test if it is == to dr in netcdf
    drw = np.empty(rw.shape) 
    drw[0] = rw[0] - 0
    drw[1:] = (rw[1:] - rw[0:-1]) * 1e6

    drd = np.empty(rd.shape) 
    drd[0] = rd[0] - 0
    drd[1:] = (rd[1:] - rd[0:-1]) * 1e6

    for t in range(data.variables['t'].shape[0]):

        g('reset')
        g('set output "' + outfolder + 'plot_spec_' + str("%03d" % t) + '.svg"')
        g('set logscale xy')
        g('set ylabel "[mg^{-1} μm^{-1}]"')
        g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
        g('set grid')
        g('set nokey')

        # FSSP range
        g('set arrow from .5,' + str(ymin) + 'to .5,' + str(ymax) + 'nohead')
        g('set arrow from 25,' + str(ymin) + 'to 25,' + str(ymax) + 'nohead')

        g('set xlabel "particle radius [μm]" ')

        nw = data.variables['wradii_m0'][t,:] / drw
        nd = data.variables['dradii_m0'][t,:] / drd
      
        plot_rw = Gnuplot.PlotItems.Data(rw, nw, with_="fsteps", title="wet radius")
        plot_rd = Gnuplot.PlotItems.Data(rd, nd, with_="fsteps", title="dry radius")

        g.plot(plot_rw, plot_rd)

def main():

    outfile = "test_spectrum.nc"
    out_bin = '{"wradii": {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0]},\
                "dradii": {"rght": 1e-6, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [0]}}'

    # run parcel run!
    parcel(dt = .5, sd_conc = 1024, outfreq = 40,  outfile = outfile, out_bin = out_bin)

    data = netcdf.netcdf_file(outfile, "r")

    # plotting 
    plot_spectrum(data, outfolder="../outputs/")

    # cleanup
    subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()

