import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf as nc
import subprocess

import parcel as pc

#checking if chemistry is indeed switched on/off 
def test_chem_off():
  outfile_on = "test_on.nc"
  pc.parcel(outfile=outfile_on, SO2_g_0=100)
  nc_on = nc.netcdf_file(outfile_on)

  outfile_off = "test_off.nc"
  pc.parcel(outfile=outfile_off, SO2_g_0=0, O3_g_0=0, H2O2_g_0=0)
  nc_off = nc.netcdf_file(outfile_off)

  assert(nc_on.SO2_g_0 > 0)
  assert(nc_off.SO2_g_0 == 0)
  assert(nc_on.variables.has_key("SO2_g"))
  assert(nc_on.variables.has_key("SO2_a"))
  assert(not nc_off.variables.has_key("SO2_g"))
  assert(not nc_off.variables.has_key("SO2_a"))

  subprocess.call(["rm", outfile_on])
  subprocess.call(["rm", outfile_off])
