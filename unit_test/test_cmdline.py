import sys,subprocess,filecmp
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
import pytest

# examples of usage, type of value should be float (??todo??)
@pytest.mark.parametrize("arg",[{"r_0" : 0.01, "T_0" : 298., "z_max" : 300.}, 
                                {"w" : 2., "p_0" : 1e5}, 
                                {"dt" : 0.5, "outfreq" : 20},
                                {"sd_conc" : 32., "kappa" : 1.}, 
                                {"SO2_g_0" : 0., "O3_g_0" : 0., "H2O2_g_0" : 0.},
                                {"gstdev" : 1.1}, 
                                {"outfreq" : 2}, 
                                {"out_bin" : ["radii:0/1/1/lin/wet/3"]},
                                {"out_bin" : ["r1:0/1/1/lin/wet/3","r2:1e-10/1e10/1/log/wet/3"]}, 
                                {"chem_dsl" : 1, "out_bin" : ["chem1:0/1/2/lin/wet/O3_a,H2O2_a", "chem1:0/1/2/lin/wet/SO2_a"]} 
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
    # arg name
    list_arg.append("--" + key)

    # arg value
    if type(value) == list:
      for v in value:
        list_arg.append(str(v))
    else:
      list_arg.append(str(value))

  subprocess.check_call(list_arg)

  # comparing if the output is the same
  subprocess.check_call(["diff", file, str(tmpdir.join("test_pyt.nc"))])
