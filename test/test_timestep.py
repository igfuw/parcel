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
import subprocess

def test_timestep(eps=0.2):
    """
    This test checks how the timestep of the simulation affects the maximum supersaturation (RH)
    and the concentration of the activated particles (N).
    The expected result is to see nearly constant RH and N for small timesteps (.001 - .03 for this setup)
    and then decrease of RH and increase of N for bigger timesteps.
    
    The testing is done by:
      - checking if the results obtained from simulations with different timesteps
        do not differ from the referential one (the one with the smallest timestep) more than eps times
      - checking if the results are the same as in the referential simulation 
        (stored in refdata folder as a gnuplot script)

    """
    outfile_nc = "test.nc"
    outfile_gp = "test_timestep_gnuplot"

    # initial values
    RH_init = .99999
    T_init  = 280.
    p_init  = 100000.
    r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))

    dt_list = [1e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 8e-3, 1e-2, 2e-2, 4e-2, 8e-2, 1e-1, 2e-1, 4e-1, 8e-1, 1.]

    # lists to store RH_max and N at the end of the simulation from each test run
    RH_list = []
    N_list  = []

    for dt in dt_list:
        print dt
        parcel(dt=dt, outfreq = 1,   outfile = outfile_nc,\
                w = 1., T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 20, \
                mean_r = 5e-8, gstdev = 1.5, n_tot = 1e9, sd_conc = 1000.)

        f_out  = netcdf.netcdf_file(outfile_nc, "r")
        RH_max = f_out.variables["RH"][:].max()
        N_end  = sum(f_out.variables["conc"][-1, -9:]) # concentration of drops > 1e-6 m
        
        RH_list.append((RH_max - 1)*100)  # [%]
        N_list.append(N_end / 1e6)        # [1/mg]
     
    # testing if the values of variables do not differ from ref. more than eps times
    for var in (RH_list, N_list):
        for val in var:
            assert abs(val - var[0]) <= eps * var[0]
 
    # save all the gnuplot commands and data needed for the plots in a file 
    Gnuplot.GnuplotOpts.prefer_inline_data = 1
    g = Gnuplot.Gnuplot(filename = outfile_gp)
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

    # do the plotting
    subprocess.call(["gnuplot", outfile_gp])   

    # check if the plots are created
    assert(os.path.isfile('plot_timestep_RH.svg'))
    assert(os.path.isfile('plot_timestep_N.svg'))
    # check if the plots stay the same
    assert(filecmp.cmp(outfile_gp, 'test/refdata/test_timestep_gnuplot_orig'))

    # cleanup
    subprocess.call(["rm", outfile_nc])
    subprocess.call(["rm", outfile_gp])
