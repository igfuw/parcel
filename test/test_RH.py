import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import pytest

#TODO: sprawdzic dla wiekszych z_max  
# wprowadzic epsilon (w sensie >= (1-eps))? 

@pytest.mark.parametrize("kappa", [0, 0.5, 0.8])
def test_RH_gt1(tmpdir, kappa):
    str_f = str(tmpdir.join("test_pcl.nc"))
    pc.parcel(outfile=str_f,  outfreq=1, kappa=kappa)
    f_out = netcdf.netcdf_file(str_f, "r")
    assert f_out.variables["RH"][-1] >= .99

#TODO: dopisac test, ktory uwzglednia krzywizne itp. w porownaniu
