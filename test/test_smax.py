import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import os
import pytest
import pdb

#checking if calculated RH_max is >= the largest value in a netcdf file  
@pytest.mark.parametrize("outfreq", [1, 2, 20])
def test_smax(tmpdir, outfreq):
    str_f = str(tmpdir.join("test_pcl.nc"))
    print str_f
    pc.parcel(outfile=str_f,  outfreq=outfreq)
    f_out = netcdf.netcdf_file(str_f, "r")
    if outfreq == 1:
        assert f_out.variables["RH"][:].max() == f_out.RH_max
    else:
        assert f_out.variables["RH"][:].max() <= f_out.RH_max


