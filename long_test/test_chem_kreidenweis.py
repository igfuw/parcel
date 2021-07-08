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
from chemical_plot import plot_chem
from kreidenweis import plot_fig1
from kreidenweis import plot_fig2
from kreidenweis import plot_fig3
from kreidenweis import plot_pH_size_dist
from thesis_profiles import thesis_profiles
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

    p_dict['sd_conc']  = 1024
    p_dict['outfreq']  = 10 / (p_dict['dt'] * p_dict['w'])

    p_dict['out_bin']  =  p_dict['out_bin'][:-1] + \
      ', "chem"  : {"rght": 1e-4, "left": 1e-6,  "drwt": "wet", "lnli": "log", "nbin": 100, "moms": ["H"]},\
         "radii" : {"rght": 1e-4, "left": 1e-6,  "drwt": "wet", "lnli": "log", "nbin": 100, "moms": [3]},\
         "acti"  : {"rght": 1e-3, "left": 1e-6,  "drwt": "wet", "lnli": "log", "nbin": 1,   "moms": [0,3]},\
         "chemd" : {"rght": 1e-6, "left": 1e-8,  "drwt": "dry", "lnli": "log", "nbin": 100, "moms": ["S_VI", "H2O2_a", "O3_a"]},\
         "specd" : {"rght": 1e-6, "left": 1e-8,  "drwt": "dry", "lnli": "log", "nbin": 100, "moms": [0, 1, 3]}}'

    pp.pprint(p_dict)

    # run parcel
    parcel(**p_dict)

    #simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")

    # removing all netcdf files after all tests                                      
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])

    request.addfinalizer(removing_files)
    return data

def test_chem_plot(data):
    """
    quicklook for chemistry
    """
    plot_chem(data, output_folder="plots/outputs", output_title='/test_chem_kreidenwies_')

def test_thesis_profiles(data):
    """
    profiles from parcel test for thesis
    """
    thesis_profiles(data, output_folder="plots/outputs")

def test_chem_fig1(data):
    """
    Fig 1 from Kreidenweis et al 2003
    """
    plot_fig1(data, output_folder="plots/outputs", output_title='/Kreidenweis_fig1')

def test_chem_fig2(data):
    """
    size distribution plot for dy radius
    """
    plot_fig2(data, output_folder="plots/outputs", output_title='/Kreidenweis_fig2')

def test_chem_fig3(data):
    """
    size distribution plot for S6
    """
    plot_fig3(data, output_folder="plots/outputs", output_title='/Kreidenweis_fig3')

def test_chem_pH_size_dist(data):
    """
    size distribution plot for S6
    """
    plot_pH_size_dist(data, output_folder="plots/outputs", output_title='/Kreidenweis_pH_todo')

def test_chem_sulfate_formation(data):
    """
    sulfate formation vs water-weighted average pH

    sulfate formed from H2O2 vs sulfate formed from O3
    """

    p = data.variables["p"][-1]
    T = data.variables["T"][-1]
    rhod = data.variables["rhod"][-1]

    # water weighted average pH at the end of the simulation
    r3     = data.variables["radii_m3"][-1,:]
    n_H    = data.variables["chem_H"][-1,:] / cm.M_H

    nom = 0

    for it, val in enumerate(r3):
        if val > 0:
            nom += (n_H[it] / (4./3 * math.pi * val * 1e3)) * val
    den  = np.sum(r3[:])                      # to liters

    pH  = -1 * np.log10(nom / den)

    print(" ")
    print(" ")
    print("water weighted average pH = ", pH, " vs 4.82-4.85 from size resolved models ")

    s6_ini = data.variables["chemd_S_VI"][0, :]
    s6_end = data.variables["chemd_S_VI"][-1, :]

    sulf_ppt = fn.mix_ratio_to_mole_frac((np.sum(s6_end) - np.sum(s6_ini)), p, cm.M_H2SO4, T, rhod) * 1e12
   
    print(" ")
    print("sulfate formation (ppt) = ", sulf_ppt, " vs 170-180 from size resolved models")

    ini_O3   = data.variables["O3_g"][0] / cm.M_O3 + \
               data.variables["chemd_O3_a"][0, :].sum() / cm.M_O3

    ini_H2O2 = data.variables["H2O2_g"][0] / cm.M_H2O2 + \
               data.variables["chemd_H2O2_a"][0, :].sum() / cm.M_H2O2


    end_O3   = data.variables["O3_g"][-1] / cm.M_O3 + \
               data.variables["chemd_O3_a"][-1, :].sum() / cm.M_O3

    end_H2O2 = data.variables["H2O2_g"][-1] / cm.M_H2O2 + \
               data.variables["chemd_H2O2_a"][-1, :].sum() / cm.M_H2O2

    sulf_ppt_H2O2 = fn.mix_ratio_to_mole_frac((np.sum(ini_H2O2) - np.sum(end_H2O2)) * cm.M_H2SO4, p, cm.M_H2SO4, T, rhod) * 1e12
    sulf_ppt_O3   = fn.mix_ratio_to_mole_frac((np.sum(ini_O3) - np.sum(end_O3)) * cm.M_H2SO4, p, cm.M_H2SO4, T, rhod) * 1e12

    print(" ")
    print("sulfate formation from H2O2 (ppt) = ", sulf_ppt_H2O2, " vs 85-105 from size resolved models")
    print("sulfate formation from O3 (ppt)   = ", sulf_ppt_O3, " vs 70-85 from size resolved models")

    print(" ")
    n_tot = data.variables["acti_m0"][12, 0] * rhod * 1e-6
    print("N of droplets           = ", n_tot, " in cm3")
    print("maximum supersaturation = ", (data.RH_max - 1) * 100, "%")
