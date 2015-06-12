import sys,subprocess,filecmp
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc

def test_cmdline(tmpdir):
  # calling from Python
  file = str(tmpdir.join("test.nc"))
  pc.parcel(outfile=file)
  # renaming the output file
  subprocess.call(["mv", file, str(tmpdir.join("test_pyt.nc"))])

  # calling via subprocess
  subprocess.call(["mv", file, str(tmpdir.join("test_pyt.nc"))])
  subprocess.check_call(["python", "parcel.py", "--outfile="+file])

  # comparing if the output is the same
  subprocess.check_call(["diff", file, str(tmpdir.join("test_pyt.nc"))])
