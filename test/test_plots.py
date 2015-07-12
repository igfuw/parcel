import sys, os
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat")
import subprocess
from scipy.io import netcdf
from parcel import parcel
from parcel_plot import RH_plot
from profiles_plot import plot_profiles
import pytest
import pdb

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


def test_RH_plot(data):
    RH_plot(data)


def test_profiles(data):
    plot_profiles(data)

