import sys
sys.path.insert(0, "../")
sys.path.insert(0, "/home/piotr/Piotr/IGF/parcel/plots/comparison")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")

from parcel import parcel
from libcloudphxx import common
from timestep_plot import timestep_plot

from scipy.io import netcdf
import numpy as np
import pytest
import os, glob
import subprocess
import pdb

"""
This set of tests checks how the timestep of the simulation affects
the maximum supersaturation (RH) and the concentration of the activated
particles (N).
The expected result is to see nearly constant RH and N for small timesteps
(.001 - .03 for this setup) and then decrease of RH and increase of N
for bigger timesteps.
"""

Dt_list = [1e-3, 2e-3, 4e-3, 8e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 1.]

# runs all simulations
# returns data with values of RH_max and N at the end of simulations
# keep the netcdf alive for all tests
@pytest.fixture(scope="module")
def data(request):
    # initial values
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
        parcel(dt=dt, outfreq = int(100/dt),   outfile = outfile_nc,\
                w = 1., T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 200, \
                sd_conc = 1000, \
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [5e-8], "gstdev": [1.5], "n_tot": [1e9]}}',
                out_bin = '{"radii": {"rght": 1, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 1e-06}}'
              )

        f_out  = netcdf.netcdf_file(outfile_nc, "r")
        RH_max = f_out.RH_max
        N_end  = f_out.variables["radii_m0"][-1,0] # concentration of drops > 1e-6 m

        RH_list.append((RH_max - 1)*100)  # [%]
        N_list.append(N_end / 1e6)        # [1/mg]

    data = {"RH" : RH_list, "N" : N_list, "dt" : Dt_list}

    # removing all netcdf files after all tests
    def removing_files():
         for file in glob.glob("timesteptest_dt*"):
            subprocess.call(["rm", file])
    request.addfinalizer(removing_files)
    return data


def test_timestep_eps(data, eps=0.01, dt_lim=0.01):
    """
    checking if the results obtained from simulations with different timesteps
    do not differ from the referential one (the one with the smallest timestep)
    more than eps times
    (Unitill we think of a better convergence test, the check is done for the
    smallest timesteps. This is done in order to avoid too big epsilon.
    """
    for var, val in data.iteritems():
        # check for RH and N
        if var in ["RH", "N"]:
            # for simulations with small timesteps
            for idx in range(len(data["dt"])):
                if data["dt"][idx] < dt_lim:
                    # assert that the results are close to the one with the smallest timestep
                    assert np.isclose(val[idx], val[0], atol=0, rtol=eps), str(val[idx]) + str(val[0])


@pytest.mark.parametrize("dt", Dt_list)
def test_timestep_diff(data, dt, eps=2.e-4):
    """
    checking if the results are close to the referential ones
    (stored in refdata folder)
    """

    filename = "timesteptest_dt=" + str(dt) + ".nc"
    f_test = netcdf.netcdf_file(filename, "r")
    f_ref  = netcdf.netcdf_file(os.path.join("/home/piotr/Piotr/IGF/parcel/long_test/refdata", filename), "r")
    for var in ["t", "z", "th_d", "T", "p", "r_v", "rhod"]:
         assert np.isclose(f_test.variables[var][:], f_ref.variables[var][:], atol=0, rtol=eps).all(), "differs e.g. " + str(var) + "; max(ref diff) = " + str(np.where(f_ref.variables[var][:] != 0., abs((f_test.variables[var][:]-f_ref.variables[var][:])/f_ref.variables[var][:]), 0.).max())


def test_timestep_plot(data):
    timestep_plot(data, output_folder="plots/outputs")
