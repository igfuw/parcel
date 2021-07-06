# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../../../")
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

# def plot_spectrum(data, data2,  outfolder):
def plot_spectrum(data, mean_r , velo, n_tot, gstdev, outfolder):

    my_path = os.path.abspath(outfolder) # Figures out the absolute path for you in case your working directory moves around.

    ymax = 1e4
    ymin = 1e-1

    rw = data.variables["wradii_r_wet"][:] * 1e6
    rd = data.variables["dradii_r_dry"][:] * 1e6

    d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10)
    d_log_rw = math.log(rw[2], 10) - math.log(rw[1], 10)

    # rw2 = data2.variables["wradii_r_wet"][:] * 1e6
    # rd2 = data2.variables["dradii_r_dry"][:] * 1e6

    #TODO - add test if it is == to dr in netcdf
    drw = np.empty(rw.shape)
    drw[0] = rw[0] - 0
    drw[1:] = (rw[1:] - rw[0:-1]) * 1e6

    drd = np.empty(rd.shape)
    drd[0] = rd[0] - 0
    drd[1:] = (rd[1:] - rd[0:-1]) * 1e6

    # drw2 = np.empty(rw2.shape)
    # drw2[0] = rw2[0] - 0
    # drw2[1:] = (rw2[1:] - rw2[0:-1]) * 1e6
    #
    # drd2 = np.empty(rd2.shape)
    # drd2[0] = rd2[0] - 0
    # drd2[1:] = (rd2[1:] - rd2[0:-1]) * 1e6

    for t in range(data.variables['t'].shape[0]):

        nd = data.variables['dradii_m0'][t,:] * data.variables["rhod"][t] / d_log_rd
        nw = data.variables['wradii_m0'][t,:] * data.variables["rhod"][t] / d_log_rw
        # nw = data.variables['wradii_m0'][t,:] / drw
        # nd = data.variables['dradii_m0'][t,:] / drd
        # nw2 = data2.variables['wradii_m0'][t,:] / drw2
        # nd2 = data2.variables['dradii_m0'][t,:] / drd2

        fig = plt.figure()
        fig.set_size_inches(18.5, 10.5)
        plt.step(rw, nw* 1e-6, label="wet radius")
        # plt.step(rw2, nw2,  label="1/10")
        plt.step(rd, nd* 1e-6,  label="dry radius")
        # plt.step(rd, nd, where='mid', label="dry radius")
        plt.xlabel("particle radius [μm]")
        # plt.xlim((0,15))
        # plt.ylim((1,160))
        plt.xscale('log')
        plt.xlabel("particle radius [μm]")
        plt.ylabel(r"$\frac{dN}{dlog_{10}(D)} [cm^{-3}]$")
        # plt.xlim((0.001,10))
        plt.legend()
        plt.grid(True,  linestyle='-.')
        fig.savefig(os.path.join(my_path, f'plot_spec_for_mean_radius={mean_r}_updraft_vel={velo}_time={t}.png'))
        plt.clf()
        plt.cla()
        plt.close()

def plot_init_spectrum(data, mean_r , velo, n_tot, gstdev,outfolder):

    """
    Plot the initial dry diameter distribution and compare it with the analitycal solution

    """
#Cumulus - TEST
    # n_tot   = 90e6
    # mean_r  = 0.24e-6
    # gstdev  = 1.58
    # n_tot2   = 15e6
    # mean_r2  = 0.24e-6
    # gstdev2  = 1.25

    my_path = os.path.abspath(outfolder) # Figures out the absolute path for you in case your working directory moves around.
    # from ncdf file attributes read out_bin parameters as a dictionary ...
    # out_bin = ast.literal_eval(getattr(data, "out_bin"))
    # ... and check if the spacing used in the test was logarithmic
    # assert out_bin["wradii"]["lnli"] == 'log', "this test should be run with logarithmic spacing of bins"

    # parcel initial condition
#    rd = data.variables["wradii_r_wet"][:] # left bin edges
    rd = data.variables["dradii_r_dry"][:] + data.variables["dradii_dr_dry"][:]/2 # left bin edges
    rw = data.variables["wradii_r_wet"][:] + data.variables["wradii_dr_wet"][:]/2


    rd2 = data.variables['dradii_m1'][0,:]
    rw2 = data.variables['wradii_m1'][0,:]
    print(rw2/rd2)
    # for comparison, model solution needs to be divided by log(d2) - log(d2)
    # since the test is run with log spacing of bins log(d2) - log(d1) = const
    # drw = np.empty(rw.shape)
    # drw[0] = rw[0] - 0
    # drw[1:] = (rw[1:] - rw[0:-1]) * 1e6
    # drd = np.empty(rd.shape)
    # drd[0] = rd[0] - 0
    # drd[1:] = (rd[1:] - rd[0:-1]) * 1e6
    d_log_rd = math.log(rd[2], 10) - math.log(rd[1], 10)
    d_log_rw = math.log(rw[2], 10) - math.log(rw[1], 10)

    # initial size distribution from the model
    model_dry = data.variables['dradii_m0'][0,:] * data.variables["rhod"][0] / d_log_rd
    model_wet = data.variables['wradii_m0'][0,:] * data.variables["rhod"][0] / d_log_rw
    # np.savetxt('TEST_for _HP.csv', [p for p in zip(rd, data.variables['wradii_m0'][0,:])], delimiter=',', fmt='%s')
    # variables for plotting theoretical solution
    def log10_distr2(mu, sigma, N, r):

        return (np.exp(-(np.log(r) - np.log(mu))**2 / (2 * np.log(sigma)**2))  / ( np.log(sigma)* np.sqrt(2 * np.pi)))*N*np.log(10)

    radii = np.logspace(-2, 2, 1000) * 1e-6
    theor = np.empty(radii.shape)
    theor2 = np.empty(radii.shape)
    for it in range(radii.shape[0]-1):
        # theor[it] = fn.log10_size_of_lnr(n_tot, mean_r, math.log(radii[it], 10), gstdev)
        theor[it] = log10_distr2(mean_r, gstdev, n_tot, radii[it])
        # theor2[it] = fn.log10_size_of_lnr(n_tot2, mean_r2, math.log(radii[it], 10), gstdev2)
    print(theor)
    fig1 = plt.figure()
    fig1.set_size_inches(18.5, 10.5)
    plt.plot(radii * 1e6,  theor  * 1e-6,  label="theory")
    # plt.plot(radii * 1e6,  theor2 * 1e-6,  label="theory2")
    plt.step(rd   * 1e6,  model_dry  * 1e-6,  label="model_dry")
    plt.step(rw    * 1e6,  model_wet  * 1e-6,  label="model_wet")
    plt.xscale('log')
    plt.xlabel("particle radius [μm]")
    plt.ylabel(r"$\frac{dN}{dlog_{10}(D)} [cm^{-3}]$")
    plt.xlim((0.001,10))
    # plt.ylim((0,400))
    plt.legend()
    plt.grid(True,  linestyle='-.')
    fig1.savefig(os.path.join(my_path, f'init_spectrum_for_mean_radius={mean_r}_updraft_vel={velo}.png'))
    # plt.savefig(os.path.join(output_folder, "mean radius= "+str(mu[j])+"Number of Super_droplets= "+str(SD_conc)+ " LWC.png"))
    plt.clf()
    plt.cla()
    plt.close()

def main():

    w_list = [1, 3, 5]
    sstp_cond = [10]
    mu = [0.24e-6, 0.06e-6, 0.54e-6, 0.14e-6]
    N_tot = [90e6, 70e6, 90e6, 80e6]
    stdev = [1.58, 1.38, 2.2, 1.2]

    # initial values
    RH_init = .98
    T_init  = 288.
    p_init  = 100000.
    r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))
    SD_conc = 1000
    H_max = 250
    col = ['red', 'blue', 'green']


    for j in range(len(mu)):
        folder_name = f"/home/piotr/Piotr/IGF/parcel3/parcel/plots/wyniki_python3/GL_II_2021/distribution_{mu[j]}"
        # X = os.rmdir(folder_name)
        folder = os.mkdir(folder_name)

        for w in w_list:

            print("updraft velosity", w, "mean value", mu[j])
            outfile = "mean radius= "+str(mu[j])+",updraft_velocity= "+ str(w)+".nc"
            # outfile2 = "test_spectrum2_TEST.nc"
            out_bin = '{"wradii": {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 100, "moms": [0,1]},\
                "dradii": {"rght": 1e-4, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 100, "moms": [0,1]}}'

            parcel(dt=1, outfreq = 1, w = w, T_0 = T_init, p_0 = p_init,\
                    r_0 = r_init, z_max = H_max, sstp_cond = 10, sd_conc = SD_conc, \
                    aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": ['+str(mu[j])+'], "gstdev": ['+str(stdev[j])+'], "n_tot": ['+str(N_tot[j])+']}}',\
                    outfile = outfile, out_bin = out_bin)
                        # math.ceil(1./w)
    # parcel(dt = 1, T_0 = T_init, p_0 = p_init, RH_0 = .98, sstp_cond =10, z_max = H, w=W,  sd_conc = SD, outfreq = 1, aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.08e-6, 0.24e-6], "gstdev": [1.78, 1.25], "n_tot": [90e6, 15e6]}}', outfile = outfile2, out_bin = out_bin)

            data = netcdf.netcdf_file(outfile, "r")
            # data2 = netcdf.netcdf_file(outfile2, "r")
    # plotting
    # plot_spectrum(data, data2, outfolder="/home/piotr/Piotr/IGF/parcel3/parcel/plots/wyniki_python3/GL_II_2021/TEST")
            # plot_init_spectrum(data, mu[j], w, N_tot[j],stdev[j],outfolder="/home/piotr/Piotr/IGF/parcel3/parcel/plots/wyniki_python3/GL_II_2021/TEST")

            plot_spectrum(data, mu[j], w, N_tot[j],stdev[j], outfolder=folder_name)


if __name__ == '__main__':
    main()
