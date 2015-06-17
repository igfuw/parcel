import sys,subprocess,filecmp
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
import pytest

# examples of usage, type of value should be float (??todo??)
@pytest.mark.parametrize("arg",[{"r_0" : 0.01, "T_0" : 298., "z_max" : 300.}, 
                                {"w" : 2., "p_0" : 1e5}, 
                                {"dt" : 1., "outfreq" : 1.},
                                {"sd_conc" : 32., "kappa" : 1.}, 
              pytest.mark.xfail({"mean_r" : .03e-6, "stdev" : 1., "n_tot" : 20.e6}),
                                {"SO2_0" : 0., "O3_0" : 0., "H2O2_0" : 0.}
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
  subprocess.check_call(["diff", file, str(tmpdir.join("test_pyt.nc"))])
