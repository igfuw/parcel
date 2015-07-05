import sys, glob, os
import subprocess
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "test/plots")
from scipy.io import netcdf
from parcel import parcel
from plot_pressure import plot_pressure_opt
import pytest
import pdb

Pprof_list = ["pprof_const_rhod", "pprof_const_th_rv",  "pprof_piecewise_const_rhod"]

@pytest.fixture(scope="module", params = [0.1, 1, 10])
def data(request):
    print "\n TWORZE data, dt = ", request.param
    data = {}
    data["dt"] = request.param
    for pprof in Pprof_list:
        filename = "test_" + pprof + str(request.param) + ".nc"
        parcel(dt=request.param, outfreq = 10, pprof = pprof, outfile=filename)
        data[pprof] = netcdf.netcdf_file(filename)

    def removing_files():
        print "\n ZABIJAM data, dt = ", request.param
        # cos mi * przy rm nie dzialala - TODO                          
        for file in glob.glob("test_pprof*"):
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
        assert abs(data[pprof].variables[var][-1] - data[pprof_ref].variables[var][-1]) <= eps * data[pprof_ref].variables[var][-1]
    # testing maximum value of RH 
    assert abs(data[pprof].variables["RH"][:].max() - data[pprof_ref].variables["RH"][:].max()) <= eps * data[pprof_ref].variables["RH"][:].max()

 
@pytest.mark.xfail #TODO                                                  
@pytest.mark.parametrize("pprof", Pprof_list)
def test_pressure_diff(data, pprof):
        filename =  "test_"+pprof+str(data["dt"])+".nc"
        filename_nc4 = filename.replace(".nc", "_nc4.nc")
        subprocess.call(["nccopy", "-k", "4", filename, filename_nc4])
        subprocess.check_call(["h5diff", "--delta=1e-18", os.path.join("test/refdata", filename_nc4), filename_nc4])

#TODO  - nie dziala tak jak chce i nie mam pojecia czemu...
#def test_pressure_plot(data):
#    try:
#    plot_pressure_opt(data)
#    except  KeyError:
        #pytest.fail("Unexpected MyError ..")
#        pass
