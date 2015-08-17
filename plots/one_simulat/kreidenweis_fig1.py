import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "../../")

from scipy.io import netcdf
import numpy as np
import pytest
import math
import subprocess

from parcel import parcel
from libcloudphxx import common as cm

def plot_fig1(data, output_folder = '', output_title = ''):

    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    spn_idx = 7
    #spn_idx = int(math.ceil(float(f_out_chem['open'].chem_spn)/float(f_out_chem['open'].outfreq)))

    # plot settings
    plt.figure(1, figsize=(16,12))
    plots = []

    for i in range(6):
      plots.append(plt.subplot(2,3,i+1))
                             #(rows, columns, number)

    plots[0].set_xlabel('lwc g/kg dry air')
    plots[0].grid()
    #plots[1].set_xlabel('SO2 conc (total-ppb)')
    plots[1].grid()
    #plots[2].set_xlabel('average pH')
    plots[2].grid()
 
    for ax in plots:
      ax.set_ylabel('t [s]')

    # read in variables
    t   = data.variables["t"][spn_idx:]
    n_H = np.sum(data.variables["chem_H"][spn_idx:], axis=1) / cm.M_H
    vol = np.sum(data.variables["radii_m3"][spn_idx:], axis=1) * 4/3. * math.pi * 1e3  #litres
    pH  = -1 * np.log10(n_H / vol)

    # calculate lwc
    plots[0].set_xticks([0., 0.5, 1, 1.5, 2, 2.5])
    plots[0].plot(np.sum(data.variables["radii_m3"][spn_idx:], axis=1) * 4. / 3 * math.pi * 998.2 * 1e3, t, "b.-")

    def aq_to_g(chem_aq, rhod, T, M, p):
         return chem_aq * rhod * cm.R * T / M  / p

    def g_to_aq(chem_g, rhod, T, M, p):
         return chem_g * M * p / rhod / cm.R / T

    nso2_g = g_to_aq(data.variables["SO2_g"][spn_idx:], data.variables["rhod"][spn_idx:],
                   data.variables["T"][spn_idx:], cm.M_SO2, data.variables["p"][spn_idx:]) / cm.M_SO2
    nso2_a = np.sum(data.variables["chem_SO2_a"][spn_idx:],  axis=1) / cm.M_SO2_H2O

    nhso3  = np.sum(data.variables["chem_HSO3_a"][spn_idx:], axis=1) / cm.M_HSO3
    nso3   = np.sum(data.variables["chem_SO3_a"][spn_idx:],  axis=1) / cm.M_SO3
    nhso4  = np.sum(data.variables["chem_HSO4_a"][spn_idx:], axis=1) / cm.M_HSO4
    nso4   = np.sum(data.variables["chem_SO4_a"][spn_idx:],  axis=1) / cm.M_SO4

    # convert aqueous phase S_IV within droplets to gase phase (mole fraction)
    S_IV_moles = nso2_g + nso2_a #+ nhso3 + nso3 + nhso4 + nso4

    tmp = 1./ data.variables["p"][spn_idx:]
    #S_IV_moles * cm.R * data.variables["rhod"][spn_idx:] / data.variables["p"][spn_idx:]

    cnv_so2_a = aq_to_g(np.sum(data.variables["chem_SO2_a"][spn_idx:],  axis=1),
                        data.variables["rhod"][spn_idx:], 
                        data.variables["T"][spn_idx:],
                        cm.M_SO2_H2O, 
                        data.variables["p"][spn_idx:])

    plots[1].set_xticks([0., 0.05, 0.1, 0.15, 0.2])
    plots[1].set_xticklabels(['0', '0.05', '0.1', '0.15', '0.2'])
    plots[1].set_xlim([0., 0.2])
    plots[1].plot(data.variables["SO2_g"][spn_idx:] * 1e9, t, "b.-")
    print (data.variables["SO2_g"][0]*1e9)
    print (data.variables["SO2_g"][1]*1e9)
    print (data.variables["SO2_g"][-1]*1e9)
    plots[1].set_xlabel('SO2 g ppb')
    plots[2].set_xticks([0., 0.002, 0.004, 0.006, 0.002])
    plots[2].plot(cnv_so2_a * 1e12, t, "b.-")
    plots[2].set_xlabel('SO2 a ppt')

    plots[3].set_xlabel('n so2 g')
    plots[3].grid()
    plots[3].plot(nso2_g, t)

    plots[4].set_xlabel('n so2 a')
    plots[4].grid()
    plots[4].plot(nso2_a, t)

    plots[5].set_xlabel('n so2 g+a')
    plots[5].grid()
    plots[5].plot(nso2_g + nso2_a, t)
 
    # TODO - calculate average pH
    #plots[2].set_xticks([3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8])
    #plots[2].plot(pH, t, "b.-")

    plt.savefig(output_folder + output_title + ".svg")
 
def main():
    # gas phase
    SO2_g_init  = 200e-12
    O3_g_init   = 50e-9
    H2O2_g_init = 500e-12
    CO2_g_init  = 360e-6
    NH3_g_init  = 100e-12
    HNO3_g_init = 100e-12

    # water vapour
    RH_init = .95
    T_init  = 285.2
    p_init  = 95000.
    r_init  = cm.eps * RH_init * cm.p_vs(T_init) / (p_init - RH_init * cm.p_vs(T_init))

    # aerosol size distribution
    mean_r = .08e-6 / 2
    gstdev = 2.
    n_tot  = 566.e6

    # process toggling
    chem_dsl = True
    chem_dsc = True
    chem_rct = True
    chem_spn = 10

    # output
    z_max       = 2400.
    dt          = .1
    w           = 1.
    outfreq     = int(z_max / dt / 100)
    sd_conc     = 1024.
    outfile     = "Kreidenweis_fig1.nc"

    # define output for moments and chemistry
    out_bin = '{"radii":   {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 500, "moms": [0, 3]},\
                "chem" :   {"rght": 1e-4, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 500,\
                               "moms": ["H", "OH", "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI"]}}'

    # run parcel, run!
#    parcel(dt = dt, z_max = z_max, w = w, outfreq = outfreq,\
#           T_0 = T_init, p_0 = p_init, r_0 = r_init,\
#           SO2_g_0 = SO2_g_init, O3_g_0 = O3_g_init, H2O2_g_0 = H2O2_g_init,\
#           CO2_g_0 = CO2_g_init, NH3_g_0 = NH3_g_init, HNO3_g_0 = HNO3_g_init,\
#           chem_sys = 'closed', outfile = outfile,\
#           sd_conc = sd_conc,\
#           chem_dsl = chem_dsl, chem_dsc = chem_dsc, chem_rct = chem_rct, chem_spn = chem_spn, \
#           out_bin = out_bin)
#
    # TODO - why do I have to repeat this import here?
    from scipy.io import netcdf
    data = netcdf.netcdf_file(outfile,   "r")

    plot_fig1(data, output_folder = "../outputs", output_title = "/Kreidenweis_fig1")

    #cleanup
    #subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()
