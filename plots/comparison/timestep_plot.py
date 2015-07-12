import sys, os, subprocess
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "../../")
from libcloudphxx import common
from scipy.io import netcdf
from parcel import parcel
import numpy as np
import Gnuplot

def timestep_plot(data, output_folder="../outputs"):

    RH_list = data['RH']
    N_list  = data['N']
    dt_list = data['dt']

    g = Gnuplot.Gnuplot()
    g('set term svg enhanced')
    g('set logscale x')
    g('set xlabel "dt [s]"')
    g('set output "' + output_folder + '/plot_timestep_RH.svg"')
    g('set ylabel "RH_{max}"')
    g.plot(Gnuplot.Data(dt_list, RH_list))
    g('set output "' + output_folder + '/plot_timestep_N.svg"')
    g('set ylabel "koncentracja koncowa [1/mg]"')
    g.plot(Gnuplot.Data(dt_list, N_list))

def main():

    Dt_list_diff = [1e-3, 2e-3, 4e-3, 8e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 1.]
    Dt_list = Dt_list_diff + [1.5e-3, 3e-3, 4e-2, 8e-2, 4e-1, 8e-1]

    # initial values
    z_max  = 200.
    w       = 1.                                        
    RH_init = .99999
    T_init  = 280.
    p_init  = 100000.
    r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))

    # lists to store RH_max and N at the end of the simulation from each test run 
    RH_list = []
    N_list  = []

    for dt in Dt_list:
        print "\nt time step", dt
        outfile_nc = "timesteptest_dt=" + str(dt) + ".nc"
        parcel(dt=dt, outfreq = int(z_max/w/dt),   outfile = outfile_nc,\
                w = w, T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = z_max, \
                mean_r = 5e-8, gstdev = 1.5, n_tot = 1e9, sd_conc = 1000., \
                radii = 1e-6 * pow(10, -3 + np.arange(26) * .2)
              )

        f_out  = netcdf.netcdf_file(outfile_nc, "r")
        RH_max = f_out.RH_max
        N_end  = sum(f_out.variables["conc"][-1, -9:]) # concentration of drops > 1e-6 m                                                                                  
        RH_list.append((RH_max - 1)*100)  # [%]                                      
        N_list.append(N_end / 1e6)        # [1/mg] 

        subprocess.call(["rm", outfile_nc])          

    data = {"RH" : RH_list, "N" : N_list, "dt" : Dt_list}
    timestep_plot(data)

if __name__ == '__main__':
    main()
