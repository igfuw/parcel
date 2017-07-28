import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import pytest

@pytest.mark.parametrize("kappa", [0.1, 0.5, 0.8])
def test_RH_gt1(tmpdir, kappa):
    """ checking if RH at the end of simulation is close to 1 
    (i.e. cloud is created) 
    """
    str_f = str(tmpdir.join("test_pcl.nc"))
    pc.parcel(outfile=str_f,  outfreq=1, aerosol='{"ammonium_sulfate": {"kappa": ' + str(kappa) + ', "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [60.0e6]}}')
    f_out = netcdf.netcdf_file(str_f, "r")
    assert f_out.variables["RH"][-1] >= .99

