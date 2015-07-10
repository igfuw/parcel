import sys, glob, os
import subprocess
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/comparison/")
from scipy.io import netcdf
from parcel import parcel
from plot_pressure import plot_pressure_opt
import numpy as np
import pytest
import pdb

Pprof_list = ["pprof_const_rhod", "pprof_const_th_rv",  "pprof_piecewise_const_rhod"]

@pytest.fixture(scope="module", params = [0.1, 1, 10])
def data(request):
    print "\n TWORZE data, dt = ", request.param
    data = {}
    data["dt"] = request.param
    for pprof in Pprof_list:
        filename = "profopttest_" + pprof + str(request.param) + ".nc"
        parcel(dt=request.param, outfreq = 10, pprof = pprof, outfile=filename)
        data[pprof] = netcdf.netcdf_file(filename)

    def removing_files():
        print "\n ZABIJAM data, dt = ", request.param
        # cos mi * przy rm nie dzialala - TODO                          
        for file in glob.glob("profopttest_pprof*"):
            subprocess.call(["rm", file])
    request.addfinalizer(removing_files)
    return data


@pytest.mark.parametrize("pprof", ["pprof_const_rhod", "pprof_const_th_rv"])
def test_pressure_opt(data, pprof, eps=0.01):
    """Testing difrent pprof options by comparing some variables to the ones 
        from reference simulation """

    # the reference option
    pprof_ref = "pprof_piecewise_const_rhod"
    # chosing variables that will be tested
    variables = ["p", "th_d", "T", "rhod", "r_v"]
    
    # testing if the values of variables do not differ from ref. more than eps times
    for var in variables:
        assert np.isclose(data[pprof].variables[var][-1], data[pprof_ref].variables[var][-1], atol=0, rtol=eps)
    # testing maximum value of RH 
    assert np.isclose(data[pprof].variables["RH"][:].max(), data[pprof_ref].variables["RH"][:].max(), atol=0, rtol=eps)

 
@pytest.mark.xfail #TODO                                                  
@pytest.mark.parametrize("pprof", Pprof_list)
def test_pressure_diff(data, pprof):
    f_ref  = netcdf.netcdf_file(os.path.join("test/refdata", 
                               "profopttest_" + pprof + str(data["dt"]) + ".nc"), "r")
    for var in f_ref.variables:
        assert np.isclose(f_ref.variables[var][:], data[pprof].variables[var][:], atol=1.e-5, rtol=0).all()


def test_pressure_plot(data):
    plot_pressure_opt(data, output_folder="plots/outputs")

