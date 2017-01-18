import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest
import copy
import pprint as pp

from parcel import parcel
from libcloudphxx import common as cm
from chem_conditions import parcel_dict
import functions as fn

@pytest.fixture(scope="module")
def data(request):

    # copy options from chem_conditions ...
    p_dict = copy.deepcopy(parcel_dict)

    # ... and modify them for the current test
    p_dict['outfile']  = "test_chem_kreidenweis.nc"

    p_dict['chem_dsl'] = True
    p_dict['chem_dsc'] = True
    p_dict['chem_rct'] = True

    p_dict['outfreq']  = 10 / (p_dict['dt'] * p_dict['w'])

    p_dict['out_bin']  =  p_dict['out_bin'][:-1] + \
      ', "chem"  : {"rght": 1e-4, "left": 1e-6,  "drwt": "wet", "lnli": "log", "nbin": 100, "moms": ["H"]},\
         "radii" : {"rght": 1e-4, "left": 1e-6,  "drwt": "wet", "lnli": "log", "nbin": 100, "moms": [3]},\
         "acti"  : {"rght": 1e-3, "left": 1e-6,  "drwt": "wet", "lnli": "log", "nbin": 1,   "moms": [0]},\
         "chemd" : {"rght": 1e-6, "left": 1e-8,  "drwt": "dry", "lnli": "log", "nbin": 100, "moms": ["S_VI" ]}}'

    # create empty lists for storing simulation results
    sd_list     = []
    pow_list    = []
    RH_max_list = []
    N_drop_list = []
    pH_list     = []
    tot_S_list  = []

    # run parcel test for super-droplet numbers between 1 and 1024
    for pw in range (11):
        sd_num = pow(2,pw)

        p_dict['sd_conc']  = sd_num
        pp.pprint(p_dict)

        # run parcel
        parcel(**p_dict)

        # save the simulation results ...
        tmp_data = netcdf.netcdf_file(p_dict['outfile'],   "r")

        # ... do the initial analysis ...
        p    = tmp_data.variables["p"][-1]
        T    = tmp_data.variables["T"][-1]
        rhod = tmp_data.variables["rhod"][-1]
        r3   = tmp_data.variables["radii_m3"][-1,:]
        n_H  = tmp_data.variables["chem_H"][-1,:] / cm.M_H
    
        # water weighted average pH at the end of the simulation:
        nom = 0
        for it, val in enumerate(r3):
            if val > 0:
                nom += (n_H[it] / (4./3 * math.pi * val * 1e3)) * val
        den  = np.sum(r3[:])                      # to liters
        pH  = -1 * np.log10(nom / den)
    
        # total sulfate formation:
        s6_ini = tmp_data.variables["chemd_S_VI"][0, :]
        s6_end = tmp_data.variables["chemd_S_VI"][-1, :]
    
        sulf_ppt = fn.mix_ratio_to_mole_frac((np.sum(s6_end) - np.sum(s6_ini)), p, cm.M_H2SO4, T, rhod) * 1e12
      
        # number of droplets at RH = max:
        n_tot = tmp_data.variables["acti_m0"][-12, 0] * rhod * 1e-6
                                        #     ^^ TODO - think of something better
        # maximum supersaturation:
        RH_max = (tmp_data.RH_max - 1) * 100

        # ... save the results ...
        sd_list.append(sd_num)
        pow_list.append(pw)
        RH_max_list.append(RH_max)
        N_drop_list.append(n_tot)
        pH_list.append(pH)
        tot_S_list.append(sulf_ppt)

        # ... and remove the simulation data files
        subprocess.call(["rm", p_dict['outfile']])

    # pack all the data into a dictionary
    data = {}
    data['number of super-droplets']  = sd_list
    data['power of two']              = pow_list
    data['cloud droplet conc. at CB'] = N_drop_list
    data['maximum supersaturation']   = RH_max_list
    data['water weighted average pH'] = pH_list
    data['total sulfate production']  = tot_S_list
 
    return data

def test_Kreidenweis_convergence(data):
    """
    test how te cloud droplet concentration, maximum supersaturation
    water weighted pH and total sulfate production depend on the 
    number of super droplets sude during computations

    """
    pp.pprint(data)

    assert np.equal(data['power of two'][-1], 10),\
       "power of two thould be 10 and is " + str(data['power of two'][-1])
    assert np.equal(data['number of super-droplets'][-1], 1024),\
       "number of super-droplets should be 1024 and is " + str(data['number of super-droplets'][-1])
    assert np.isclose(data['cloud droplet conc. at CB'][-1], 273 , atol=2),\
       " cloud droplet conc. at CB should be 273 and is " + str(data['cloud droplet conc. at CB'][-1])
    assert np.isclose(data['maximum supersaturation'][-1], 0.27, atol=0.01),\
       "maximum supersaturation should be 0.27 and is " + str(data['maximum supersaturation'][-1])
    assert np.isclose(data['water weighted average pH'][-1], 4.83, atol=0.01),\
       "water weighted average pH should be 4.83 and is " + str(data['water weighted average pH'][-1])
    assert np.isclose(data['total sulfate production'][-1], 170, atol=1),\
       "total sulfate production should be " + str(data['total sulfate production'][-1])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties

    # plot the results
    plt.figure(1)
    plt.figure(figsize=(25,20))
    plt.rcParams.update({'font.size': 30})

    plots = []
    for i in range(4):
      plots.append(plt.subplot(2,2,i+1))
      plots[i].grid()
      plots[i].set_xlabel('$\mathrm{ln_2(number \, of \, super-droplets)}$')

    for ax in plots:
      #ax.set_ylim([0, 2400])
      #ax.set_yticks([0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400])
      ax.tick_params(axis='x', pad=15)
      ax.tick_params(axis='y', pad=15)

    plots[0].plot(data['power of two'], data['cloud droplet conc. at CB'], "b.-", ms=15, lw=3.)
    plots[1].plot(data['power of two'], data['maximum supersaturation'],   "b.-", ms=15, lw=3.)
    plots[2].plot(data['power of two'], data['water weighted average pH'], "b.-", ms=15, lw=3.)
    plots[3].plot(data['power of two'], data['total sulfate production'],  "b.-", ms=15, lw=3.)

    plots[0].set_ylabel('droplet conc. at cloud base $\mathrm{[cm^{-3}]}$')
    plots[1].set_ylabel('max. supersat. [%]')
    plots[2].set_ylabel('average pH')
    plots[3].set_ylabel('total sulfate production [ppt]')

    plt.savefig("plots/outputs/Kreidenweis_convergence.pdf")
