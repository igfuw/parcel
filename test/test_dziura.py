import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
from scipy.io import netcdf
from parcel import parcel
import pytest
import Gnuplot

def test_dziura():

    RH_init = .99999
    T_init  = 280.
    p_init  = 100000.
    r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))

    f_res = open('test_dziura.txt', 'w')

    dt_list = [8e-4, 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 4e-3, 8e-3, 1e-2, 2e-2, 4e-2, 8e-2, 1e-1, 2e-1, 4e-1, 8e-1, 1.]

    for dt in dt_list:
        print dt
        parcel(dt=dt, outfreq = 1,   outfile="test.nc",\
                w = 1., T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 20, \
                mean_r = 5e-8, gstdev = 1.5, n_tot = 1e9, sd_conc = 1000., kappa = 1. )

        f_out = netcdf.netcdf_file("test.nc", "r")
        RH_max = f_out.variables["RH"][:].max()
        N_end = f_out.variables["conc"][-1, -9:]

        f_res.write(str(dt) + " " + str(RH_max) + " " + str(sum(N_end)) + "\n")
       
    g = Gnuplot.Gnuplot()
    g('set term svg enhanced')

#    g('set output "tmp_dziura.svg"')
#    z = f_out.variables['z'][:]
#    rh = (f_out.variables['RH'][:] - 1) * 100
#    g.plot(Gnuplot.Data(rh, z))

    g('set logscale x')
    g('set xlabel "dt [s]"')

    g('set output "plot_dziura1.svg" ')
    g('set ylabel "RH_{max}"')
    g('plot "test_dziura.txt" u ($1):(($2)-1)*100  notitle')

    g('set output "plot_dziura2.svg" ')
    g('set ylabel "koncentracja koncowa [1/mg]"')
    g('plot "test_dziura.txt" u ($1):($3)/1e6 notitle')
