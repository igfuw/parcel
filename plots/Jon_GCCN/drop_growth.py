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


def plot_growth_GCCN(data, output_folder = '', output_title = ''):
    import Gnuplot

    z = data.variables["z"][:]
    curr_w = data.variables["curr_w"][:]
    rd1um_mean_r =  data.variables["rd1.0um_wet_mom1"][:,0]   / data.variables["rd1.0um_wet_mom0"][:,0]   
    rd2um_mean_r =  data.variables["rd2.0um_wet_mom1"][:,0]   / data.variables["rd2.0um_wet_mom0"][:,0]   
    rd3um_mean_r =  data.variables["rd3.0um_wet_mom1"][:,0]   / data.variables["rd3.0um_wet_mom0"][:,0]   
    rd4um_mean_r =  data.variables["rd4.0um_wet_mom1"][:,0]   / data.variables["rd4.0um_wet_mom0"][:,0]   
    rd5um_mean_r =  data.variables["rd5.0um_wet_mom1"][:,0]   / data.variables["rd5.0um_wet_mom0"][:,0]   
    rd6um_mean_r =  data.variables["rd6.0um_wet_mom1"][:,0]   / data.variables["rd6.0um_wet_mom0"][:,0]   
    rd7um_mean_r =  data.variables["rd7.0um_wet_mom1"][:,0]   / data.variables["rd7.0um_wet_mom0"][:,0]   
#    rd8um_mean_r =  data.variables["rd8.0um_wet_mom1"][:,0]   / data.variables["rd8.0um_wet_mom0"][:,0]   
    rd9um_mean_r =  data.variables["rd9.0um_wet_mom1"][:,0]   / data.variables["rd9.0um_wet_mom0"][:,0]   

    g = Gnuplot.Gnuplot()# persist=1)
    g('set terminal svg')
    g('set samples 500')
    g('set key bottom')
    g('set nokey')

    g('set output "' + output_folder + output_title + '.svg"')
    g('set xrange [0:100]')
    g('set yrange [-300:400]')

    
    plot_rd1um_ud = Gnuplot.PlotItems.Data(rd1um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=1um")
    plot_rd2um_ud = Gnuplot.PlotItems.Data(rd2um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=2um")
    plot_rd3um_ud = Gnuplot.PlotItems.Data(rd3um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=3um")
    plot_rd4um_ud = Gnuplot.PlotItems.Data(rd4um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=4um")
    plot_rd5um_ud = Gnuplot.PlotItems.Data(rd5um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=5um")
    plot_rd6um_ud = Gnuplot.PlotItems.Data(rd6um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=6um")
    plot_rd7um_ud = Gnuplot.PlotItems.Data(rd7um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=7um")
#    plot_rd8um_ud = Gnuplot.PlotItems.Data(rd8um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=8um")
    plot_rd9um_ud = Gnuplot.PlotItems.Data(rd9um_mean_r  [curr_w>0] * 1e6, z[curr_w>0]-300., with_="l lc 1", title="rd=9um")

    plot_rd1um_dd = Gnuplot.PlotItems.Data(rd1um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=1um")
    plot_rd2um_dd = Gnuplot.PlotItems.Data(rd2um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=2um")
    plot_rd3um_dd = Gnuplot.PlotItems.Data(rd3um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=3um")
    plot_rd4um_dd = Gnuplot.PlotItems.Data(rd4um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=4um")
    plot_rd5um_dd = Gnuplot.PlotItems.Data(rd5um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=5um")
    plot_rd6um_dd = Gnuplot.PlotItems.Data(rd6um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=6um")
    plot_rd7um_dd = Gnuplot.PlotItems.Data(rd7um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=7um")
#    plot_rd8um_dd = Gnuplot.PlotItems.Data(rd8um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=8um")
    plot_rd9um_dd = Gnuplot.PlotItems.Data(rd9um_mean_r  [curr_w<0] * 1e6, z[curr_w<0]-300., with_="l lc 2", title="rd=9um")

    g.plot(plot_rd1um_ud, \
           plot_rd2um_ud, \
           plot_rd3um_ud, \
           plot_rd4um_ud, \
           plot_rd5um_ud, \
           plot_rd6um_ud, \
           plot_rd7um_ud, \
#           plot_rd8um_ud, \
           plot_rd9um_ud, \
           plot_rd1um_dd, \
           plot_rd2um_dd, \
           plot_rd3um_dd, \
           plot_rd4um_dd, \
           plot_rd5um_dd, \
           plot_rd6um_dd, \
           plot_rd7um_dd, \
#           plot_rd8um_dd, \
           plot_rd9um_dd
     )

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
    g('set terminal svg')
    g('set samples 500')
    g('set key bottom')
    g('set nokey')

    g('set output "' + output_folder + output_title + '.svg"')
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
    rd7um_mean_r =  data.variables["rd7.0um_wet_mom1"][:,0]   / data.variables["rd7.0um_wet_mom0"][:,0]   
    rd9um_mean_r =  data.variables["rd9.0um_wet_mom1"][:,0]   / data.variables["rd9.0um_wet_mom0"][:,0]   

    g = Gnuplot.Gnuplot()# persist=1)
    g('set output "' + output_folder + output_title + '.svg"')
    g('set terminal svg')
    g('set samples 500')
    g('set multiplot layout 3, 1')
#    g('set key bottom')


    
    plot_S = Gnuplot.PlotItems.Data(S, z-300., with_="l", title="S[%]")
    plot_S_libcloud = Gnuplot.PlotItems.Data(S_libcloud, z-300., with_="l", title="S libcloud[%]")
    plot_St = Gnuplot.PlotItems.Data(t, S, with_="l", title="S[%]")
    plot_S_libcloudt = Gnuplot.PlotItems.Data(t, S_libcloud, with_="l", title="S libcloud[%]")
    plot_zt = Gnuplot.PlotItems.Data(z-300., t, with_="l", title="z[m]")
    plot_r_vt = Gnuplot.PlotItems.Data(r_v, t, with_="l", title="rv[1]")
    plot_rhodt = Gnuplot.PlotItems.Data(rhod, t, with_="l", title="rhod")
    plot_pt = Gnuplot.PlotItems.Data(p, t, with_="l", title="p")
    plot_rd7t = Gnuplot.PlotItems.Data(t, rd7um_mean_r * 1e6,  with_="l", title="rd7um")
    plot_rd9t = Gnuplot.PlotItems.Data(t, rd9um_mean_r * 1e6,  with_="l", title="rd9um")

#    g('set xrange [-1:1]')
#    g('set yrange [-350:350]')
#    g.plot(plot_S, plot_S_libcloud)
#    g('set yrange [0:*]')
    g('set yrange [-0.05:-0.045]')
    g('set xrange [1650:1750]')
    g.plot(plot_St, plot_S_libcloudt)
#    g('set xrange [-350:350]')
#    g.plot(plot_zt)
#    g('set xrange [*:*]')
#    g.plot(plot_r_vt)
#    g.plot(plot_rhodt)
#    g.plot(plot_pt)
    g('set yrange [*:*]')
    g.plot(plot_rd9t)
    g.plot(plot_rd7t)



def main():

    data = netcdf.netcdf_file('Jon_GCCN.nc',   "r")

    # do the plotting
    plot_growth_noGCCN(data, output_folder = "../outputs", output_title = "/Jon_GCCN_drop_growth_noGCCN")
    plot_growth_GCCN(data, output_folder = "../outputs", output_title = "/Jon_GCCN_drop_growth_GCCN")
    plot_supersat(data, output_folder = "../outputs", output_title = "/Jon_GCCN_supersat")

    # cleanup
#    subprocess.call(["rm", p_dict['outfile']])

if __name__ == '__main__':
    main()
