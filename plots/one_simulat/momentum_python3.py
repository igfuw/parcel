# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../../")
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")

import ast
import functions as fn
from parcel import parcel
from libcloudphxx import common
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.io import netcdf
import numpy as np
import pytest
import os, glob
import subprocess
import pdb
import math
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.animation import FuncAnimation

def plot_spectrum(data, outfolder):


    my_path = os.path.abspath(outfolder)
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

        fig = plt.figure()
        fig.set_size_inches(18.5, 10.5)
        plt.step(rw, nw, label="1")
        plt.step(rw2, nw2,  label="1/10")
        plt.step(rd2, nd2,  label="dry radius")
        # plt.step(rd, nd, where='mid', label="dry radius")
        plt.xlabel("particle radius [μm]")
        plt.xlim((0,15))
        plt.ylim((1,160))
        plt.legend()
        plt.grid(True,  linestyle='-.')
        fig.savefig(os.path.join(my_path, 'plot_spec_' + str("%03d" % t) + '_DYCOMS.svg'))
        plt.clf()
        plt.cla()
        plt.close()

def main():

    outfile = "test_spectrum.nc"
    out_bin = '"cloud": {"rght": 2.5e-05, "moms": [0, 1, 2, 3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-07}'

    # run parcel run!
    parcel(dt = .5, sd_conc = 1024, outfreq = 40,  outfile = outfile, out_bin = out_bin)

    data = netcdf.netcdf_file(outfile, "r")

    # plotting
    plot_spectrum(data, outfolder="/home/piotr/Piotr/IGF/parcel3/parcel/wyniki_laptop")

    # cleanup
    subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()
