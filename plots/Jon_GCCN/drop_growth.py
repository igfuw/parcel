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


def plot_growth_noGCCN(data, output_folder = '', output_title = ''):
    import Gnuplot

    z = data.variables["z"][:]
    curr_w = data.variables["curr_w"][:]
    rd20_mean_r =  data.variables["rd20nm_wet_mom1"][:,0]   / data.variables["rd20nm_wet_mom0"][:,0]   
    rd31_mean_r =  data.variables["rd31nm_wet_mom1"][:,0]   / data.variables["rd31nm_wet_mom0"][:,0]
    rd152_mean_r = data.variables["rd152nm_wet_mom1"][:,0]  / data.variables["rd152nm_wet_mom0"][:,0]
    rd337_mean_r = data.variables["rd337nm_wet_mom1"][:,0]  / data.variables["rd337nm_wet_mom0"][:,0]
    rd500_mean_r = data.variables["rd500nm_wet_mom1"][:,0]  / data.variables["rd500nm_wet_mom0"][:,0]

    g = Gnuplot.Gnuplot()# persist=1)
    g('set terminal pdf')
    g('set samples 500')
    g('set key bottom')
    g('set nokey')

    g('set output "' + output_folder + output_title + '.pdf"')
    g('set xrange [0:15]')
    g('set yrange [-300:400]')

    
    plot_rd20_ud = Gnuplot.PlotItems.Data(rd20_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=20nm")
    plot_rd31_ud = Gnuplot.PlotItems.Data(rd31_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=31nm")
    plot_rd152_ud = Gnuplot.PlotItems.Data(rd152_mean_r[curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=152nm")
    plot_rd337_ud = Gnuplot.PlotItems.Data(rd337_mean_r[curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=337nm")
    plot_rd500_ud = Gnuplot.PlotItems.Data(rd500_mean_r[curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=500nm")

    plot_rd20_dd = Gnuplot.PlotItems.Data(rd20_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=20nm")
    plot_rd31_dd = Gnuplot.PlotItems.Data(rd31_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=31nm")
    plot_rd152_dd = Gnuplot.PlotItems.Data(rd152_mean_r[curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=152nm")
    plot_rd337_dd = Gnuplot.PlotItems.Data(rd337_mean_r[curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=337nm")
    plot_rd500_dd = Gnuplot.PlotItems.Data(rd500_mean_r[curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=500nm")

    g.plot(plot_rd20_ud, plot_rd31_ud, plot_rd152_ud, plot_rd337_ud, plot_rd500_ud \
          ,plot_rd20_dd, plot_rd31_dd, plot_rd152_dd, plot_rd337_dd, plot_rd500_dd \
     )

def plot_supersat(data, output_folder = '', output_title = ''):
    import Gnuplot

    z = data.variables["z"][:]
    t = data.variables["t"][:]
    S = (data.variables["RH"][:] - 1.) * 100
    r_v = data.variables["r_v"][:]
    p = data.variables["p"][:]
    rhod = data.variables["rhod"][:]
    S_libcloud = (data.variables["RH_libcloud"][:] - 1.) * 100

    g = Gnuplot.Gnuplot()# persist=1)
    g('set output "' + output_folder + output_title + '.pdf"')
    g('set terminal pdf')
    g('set samples 500')
    g('set multiplot layout 2, 3')
#    g('set key bottom')


    
    plot_S = Gnuplot.PlotItems.Data(S, z-300., with_="l", title="S[%]")
    plot_S_libcloud = Gnuplot.PlotItems.Data(S_libcloud, z-300., with_="l", title="S libcloud[%]")
    plot_St = Gnuplot.PlotItems.Data(S, t, with_="l", title="S[%]")
    plot_zt = Gnuplot.PlotItems.Data(z-300., t, with_="l", title="z[m]")
    plot_r_vt = Gnuplot.PlotItems.Data(r_v, t, with_="l", title="rv[1]")
    plot_rhodt = Gnuplot.PlotItems.Data(rhod, t, with_="l", title="rhod")
    plot_pt = Gnuplot.PlotItems.Data(p, t, with_="l", title="p")

    g('set xrange [-1:1]')
    g('set yrange [-350:350]')
    g.plot(plot_S, plot_S_libcloud)
    g('set yrange [0:*]')
    g.plot(plot_St)
    g('set xrange [-350:350]')
    g.plot(plot_zt)
    g('set xrange [*:*]')
    g.plot(plot_r_vt)
    g.plot(plot_rhodt)
    g.plot(plot_pt)



def main():

    data = netcdf.netcdf_file('Jon_GCCN.nc',   "r")

    # do the plotting
    plot_growth_noGCCN(data, output_folder = "../outputs", output_title = "/Jon_GCCN_drop_growth_noGCCN")
    plot_supersat(data, output_folder = "../outputs", output_title = "/Jon_GCCN_supersat")

    # cleanup
#    subprocess.call(["rm", p_dict['outfile']])

if __name__ == '__main__':
    main()
