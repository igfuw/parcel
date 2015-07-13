import sys, glob, os
import subprocess
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")
from scipy.io import netcdf
import numpy as np
import pytest
import pdb

from parcel import parcel
from pressure_plot import plot_pressure_opt

"""
 This set of test checks how the option pprof (choosing the method of calculating
pressure  profiles) affects results. 
"""

Pprof_list = ["pprof_const_rhod", "pprof_const_th_rv",  "pprof_piecewise_const_rhod"]

# runs all simulations 
# returns opened netcdfile files
@pytest.fixture(scope="module", params = [0.1])
def data(request):
    data = {}
    data["dt"] = request.param
    for pprof in Pprof_list:
        filename = "profopttest_" + pprof + str(request.param) + ".nc"
        parcel(dt=request.param, outfreq = 100, pprof = pprof, outfile=filename)
        data[pprof] = netcdf.netcdf_file(filename)
    # removing all netcdf files after all tests
    def removing_files():
        for file in glob.glob("profopttest_pprof*"):
            subprocess.call(["rm", file])
    request.addfinalizer(removing_files)
    return data


@pytest.mark.parametrize("pprof", ["pprof_const_rhod", "pprof_const_th_rv"])
def test_pressure_opt(data, pprof, eps=0.01):
    """    
    checking if the results obtained from simulations with different pprof      
    do not differ from the referential one  more than eps times
    """
    # the reference option
    pprof_ref = "pprof_piecewise_const_rhod"
    # chosing variables that will be tested
    variables = ["p", "th_d", "T", "rhod", "r_v"]
    
    # testing if the values of variables do not differ from ref. more than eps times
    for var in variables:
        assert np.isclose(data[pprof].variables[var][-1], data[pprof_ref].variables[var][-1], atol=0, rtol=eps)
    # testing maximum value of RH 
    assert np.isclose(data[pprof].variables["RH"][:].max(), data[pprof_ref].variables["RH"][:].max(), atol=0, rtol=eps)

 
@pytest.mark.parametrize("pprof", Pprof_list)
def test_pressure_diff(data, pprof, eps=1.e-8):
    """     
    checking if the results for all pprof option are close to the referential ones
    (stored in refdata folder)                                             
    """

    f_ref  = netcdf.netcdf_file(os.path.join("test/refdata", 
                             "profopttest_" + pprof + str(data["dt"]) + ".nc"), "r")
    for var in ["t", "z", "th_d", "T", "p", "r_v", "rhod", "RH"]:
        assert np.isclose(f_ref.variables[var][:], data[pprof].variables[var][:], atol=0, rtol=eps).all(), "differs e.g. " + str(var) + "; max(ref diff) = " + str(np.where(f_ref.variables[var][:] != 0., abs((data[pprof].variables[var][:]-f_ref.variables[var][:])/f_ref.variables[var][:]), 0.).max())


def test_pressure_plot(data):
    """
    checking if the plot function works correctly
    returning the plot showing differences between simulations 
    with different pprof options
    """
    plot_pressure_opt(data, output_folder="plots/outputs")

