import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import os
import pytest
import pdb

def test_smax(tmpdir):
    str_f = str(tmpdir.join("test_pcl.nc"))
    print str_f
    pc.parcel(outfile=str_f,  outfreq=1)
    f_out = netcdf.netcdf_file(str_f, "r")
    assert f_out.variables["RH"][:].max() == f_out.RH_max



