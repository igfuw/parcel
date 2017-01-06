import sys, os
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat")
import subprocess
from scipy.io import netcdf
import pytest

from parcel import parcel
from thesis_profiles import thesis_profiles
"""
sets of test checking if ploting functions work correctly
"""

@pytest.fixture(scope="module")
def data(request):
    parcel(outfile="file_test.nc")
    fnc = netcdf.netcdf_file("file_test.nc")
    def removing_files():
        fnc.close()
        subprocess.call(["rm", "file_test.nc"])
    request.addfinalizer(removing_files)
    return fnc

def test_thesis_profiles(data):
    thesis_profiles(data, output_folder="plots/outputs")

