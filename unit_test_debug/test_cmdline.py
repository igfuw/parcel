import sys,subprocess,filecmp
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
import pytest

# examples of usage, type of value should be float (??todo??)
@pytest.mark.parametrize("arg",[{"r_0" : 0.01, "T_0" : 298., "z_max" : 300.}, 
                                {"w" : 2., "p_0" : 1e5}, 
                                {"dt" : 0.5, "outfreq" : 20},
                                {"sd_conc" : 32}, 
                                {"SO2_g" : 0., "O3_g" : 0., "H2O2_g" : 0.},
                                {"aerosol": '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.02e-6, 0.07e-7], "gstdev": [1.4, 1.2], "n_tot": [120.0e6, 80.0e6]},' 
                                             '"gccn"            : {"kappa": 1.28, "mean_r": [2e-6], "gstdev": [1.6], "n_tot": [1e3]}}'},
                                {"outfreq" : 2},
                                {"out_bin" : '{"radii": {"rght": 1, "moms": [3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0}}'},
                                {"out_bin" : '{"r1": {"rght": 1, "moms": [3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0}, "r2": {"rght": 10000000000.0, "moms": [3], "drwt": "wet", "nbin": 1, "lnli": "log", "left": 1e-10}}'}, 
                                {"chem_dsl" : 1, "out_bin" : '{"chem1": {"rght": 1, "moms": ["O3_a", "H2O2_a"], "drwt": "wet", "nbin": 2, "lnli": "lin", "left": 0}, "chem2": {"rght": 1, "moms": ["SO2_a"], "drwt": "wet", "nbin": 2, "lnli": "lin", "left": 0}}'} 
                                ])
def test_cmdline(tmpdir, arg):
  # calling from Python
  file = str(tmpdir.join("test.nc"))
  pc.parcel(outfile=file, **arg)
  # renaming the output file
  subprocess.call(["mv", file, str(tmpdir.join("test_pyt.nc"))])

  # calling via subprocess
  # creating a list with all provided arguments
  list_arg = ["python", "parcel.py", "--outfile="+file]
  for key, value in arg.items():
    list_arg.append("--" + key + "=" + str(value))

  subprocess.check_call(list_arg)

  # comparing if the output is the same
  # this test might fail if libcloudph++ was compiled with -Ofast flag (should work with -O3)
  subprocess.check_call(["diff", file, str(tmpdir.join("test_pyt.nc"))])
