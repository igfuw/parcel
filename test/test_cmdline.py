import sys,subprocess,filecmp
sys.path.insert(0, "../")
import parcel

def test_cmdline(tmpdir):
  # calling from Python
  file_python = str(tmpdir.join("python.nc"))
  parcel.parcel(outfile=file_python)

  # calling via subprocess
  file_sbproc = str(tmpdir.join("sbproc.nc"))
  subprocess.check_call(["parcel.py", "--outfile="+file_sbproc])

  # comparing if the output is the same
  # TODO!
