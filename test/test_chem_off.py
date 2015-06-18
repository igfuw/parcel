import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf as nc

def test_chem_off(tmpdir):
  file = str(tmpdir.join("test_on.nc"))
  pc.parcel(outfile=file, SO2_0=100)
  nc_on = nc.netcdf_file(file)

  file = str(tmpdir.join("test_off.nc"))
  pc.parcel(outfile=file, SO2_0=0, O3_0=0, H2O2_0=0)
  nc_off = nc.netcdf_file(file)

  assert(nc_on.SO2_0 > 0)
  assert(nc_off.SO2_0 == 0)
  assert(nc_on.variables.has_key("SO2"))
  assert(not nc_off.variables.has_key("SO2"))
