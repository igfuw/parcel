import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf
from parcel import parcel
import pytest



@pytest.mark.parametrize("dt", [1, 10])
@pytest.mark.parametrize("pprof", ["pprof_const_rhod", "pprof_const_th_rv"])
def test_pressure_opt(dt, pprof, eps=0.01):
    """Testing difrent pprof options by comparing some variables to the ones 
        from reference simulation """

    # the reference option
    pprof_ref = "pprof_piecewise_const_rhod"
    # running parcel model for the reference option and testing one
    parcel(dt=dt, outfreq = 10, pprof = pprof, outfile="test_" + pprof + ".nc")
    parcel(dt=dt, outfreq = 10, pprof = pprof_ref, outfile="test_" + pprof_ref + ".nc")
    # opening netcdf files
    f_test = netcdf.netcdf_file("test_"+pprof+".nc", "r")
    f_ref = netcdf.netcdf_file("test_"+pprof_ref+".nc", "r")

    # chosing variables that will be tested
    variables = ["p", "th_d", "T", "rhod", "r_v"]
    
    # testing if the values of variables do not differ from ref. more than eps times
    for var in variables:
        assert abs(f_test.variables[var][-1] - f_ref.variables[var][-1]) <= eps * f_ref.variables[var][-1]
    # testing maximum value of RH 
    assert abs(f_test.variables["RH"][:].max() - f_ref.variables["RH"][:].max()) <= eps * f_ref.variables["RH"][:].max()

 
