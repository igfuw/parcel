import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from scipy.io import netcdf as nc
import subprocess

import parcel as pc

#checking if chemistry is indeed switched on/off 
def test_chem_off():
  outfile_on = "test_on.nc"
  pc.parcel(outfile=outfile_on, chem_dsl = True, chem_dsc = False, chem_rct = False,
            SO2_g = 200e-12, O3_g = 50e-9, H2O2_g = 500e-12, 
            CO2_g = 360e-6,  NH3_g = 100e-12, HNO3_g = 100e-12
           )
  nc_on = nc.netcdf_file(outfile_on)

  outfile_off = "test_off.nc"
  pc.parcel(outfile=outfile_off, chem_dsl = False, chem_dsc = False, chem_rct = False)
  nc_off = nc.netcdf_file(outfile_off)

  assert(nc_on.SO2_g > 0)
  assert(nc_off.SO2_g == 0)
  assert(nc_on.variables.has_key("SO2_g"))
  assert(nc_on.variables.has_key("SO2_a"))
  assert(not nc_off.variables.has_key("SO2_g"))
  assert(not nc_off.variables.has_key("SO2_a"))

  subprocess.call(["rm", outfile_on])
  subprocess.call(["rm", outfile_off])
