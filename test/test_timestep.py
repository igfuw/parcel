import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
from scipy.io import netcdf
from parcel import parcel
import pytest
import Gnuplot
import os
import filecmp

def test_timestep(eps=0.2):

    # initail values
    RH_init = .99999
    T_init  = 280.
    p_init  = 100000.
    r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))

    dt_list = [7e-4, 8e-4, 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 4e-3, 8e-3, 1e-2, 2e-2, 4e-2, 8e-2, 1e-1, 2e-1, 4e-1, 8e-1, 1.]
    #dt_list = [1e-1, 2e-1, 4e-1, 8e-1, 1.]

    # lists to store RH_max and N at the end of the simulation from each test run
    RH_list = []
    N_list  = []

    for dt in dt_list:
        print dt
        parcel(dt=dt, outfreq = 1,   outfile="test.nc",\
                w = 1., T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 20, \
                mean_r = 5e-8, gstdev = 1.5, n_tot = 1e9, sd_conc = 1000.)

        f_out  = netcdf.netcdf_file("test.nc", "r")
        RH_max = f_out.variables["RH"][:].max()
        N_end  = sum(f_out.variables["conc"][-1, -9:]) # concentration of drops > 1e-6 m
        
        RH_list.append((RH_max - 1)*100)
        N_list.append(N_end / 1e6)
     
    # testing if the values of variables do not differ from ref. more than eps times
    for var in (RH_list, N_list):
        for val in var:
            assert abs(val - var[0]) <= eps * var[0]
  
    g = Gnuplot.Gnuplot()
    g('set term svg enhanced')
    g('set logscale x')
    g('set xlabel "dt [s]"')

    g('set output "plot_timestep_RH.svg" ')
    g('set ylabel "RH_{max}"')
    g.plot(Gnuplot.Data(dt_list, RH_list))

    g('set output "plot_timestep_N.svg"')
    g('set ylabel "koncentracja koncowa [1/mg]"')
    g.plot(Gnuplot.Data(dt_list, N_list))

    del g # necessary to create the svg files before comparison

    assert(os.path.isfile('plot_timestep_RH.svg'))
    assert(os.path.isfile('plot_timestep_N.svg'))
    # checking if the plots are the same
    assert(filecmp.cmp('plot_timestep_RH.svg', 'test/refdata/plot_timestep_RH_orig.svg'))
    assert(filecmp.cmp('plot_timestep_N.svg',  'test/refdata/plot_timestep_N_orig.svg'))
