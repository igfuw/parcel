import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import pytest

@pytest.mark.parametrize("kappa", [0.1, 0.5, 0.8])
def test_RH_gt1(tmpdir, kappa):
    str_f = str(tmpdir.join("test_pcl.nc"))
    pc.parcel(outfile=str_f,  outfreq=1, kappa=kappa)
    f_out = netcdf.netcdf_file(str_f, "r")
    assert f_out.variables["RH"][-1] >= .99

