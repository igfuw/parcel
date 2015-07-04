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

@pytest.fixture(scope="module", params = [0.1, 1, 10])
def data_all_pressure_opt(request):
    print "\n TWORZE data, dt = ", request.param
    data = {}
    data["dt"] = request.param
    for pprof in ["pprof_const_rhod", "pprof_const_th_rv", "pprof_piecewise_const_rhod"]:
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
def test_pressure_opt(data_all_pressure_opt, pprof, eps=0.01):
    """Testing difrent pprof options by comparing some variables to the ones 
        from reference simulation """

    #pdb.set_trace()
    # the reference option
    pprof_ref = "pprof_piecewise_const_rhod"
    # opening netcdf files
    f_test = data_all_pressure_opt[pprof]
    f_ref = data_all_pressure_opt[pprof_ref]

    # chosing variables that will be tested
    variables = ["p", "th_d", "T", "rhod", "r_v"]
    
    # testing if the values of variables do not differ from ref. more than eps times
    for var in variables:
        assert abs(f_test.variables[var][-1] - f_ref.variables[var][-1]) <= eps * f_ref.variables[var][-1]
    # testing maximum value of RH 
    assert abs(f_test.variables["RH"][:].max() - f_ref.variables["RH"][:].max()) <= eps * f_ref.variables["RH"][:].max()

 
#@pytest.mark.xfail #TODO                                                  
@pytest.mark.parametrize("pprof", ["pprof_const_rhod", "pprof_const_th_rv", "pprof_piecewise_const_rhod"])
def test_pressure_diff(data_all_pressure_opt, pprof):
        filename =  "test_"+pprof+str(data_all_pressure_opt["dt"])+".nc"
        filename_nc4 = filename.replace(".nc", "_nc4.nc")
        subprocess.call(["nccopy", "-k", "4", filename, filename_nc4])
        subprocess.check_call(["h5diff", "--delta=1e-18", os.path.join("test/refdata", filename_nc4), filename_nc4])


#def test_pressure_plot(data_all_pressure_opt):
#    try:
#    plot_pressure_opt(data_all_pressure_opt)
#    except  KeyError:
        #pytest.fail("Unexpected MyError ..")
#        pass
