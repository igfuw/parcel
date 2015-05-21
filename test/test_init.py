import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import pytest

def test_supersat(tmpdir):
    str_f = str(tmpdir.join("test_pcl.nc"))
    with pytest.raises(Exception) as excinfo:
        pc.parcel(outfile=str_f,  r_0=.1)

