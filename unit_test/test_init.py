import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import pytest

"""
set of test checking arguments and initial conditions
"""

def test_supersat(tmpdir):
    """ checking if parcel rises an exception when RH_init>0 """
    str_f = str(tmpdir.join("test_pcl.nc"))
    with pytest.raises(Exception) as excinfo:
        pc.parcel(outfile=str_f,  r_0=.1)


@pytest.mark.parametrize("arg",[{"T_0" : 255},
                                {"RH_0": 1.00001},
                                {"w" : -1},
                                {"aerosol" : '{"zero_kappa": {"kappa" : 0., "mean_r": [2e-8], "gstdev": [1.2], "n_tot": [60e6]}}'},
                                {"aerosol" : '{"unity_gstdev": {"kappa" : 0.61, "mean_r": [2e-8], "gstdev": [1.], "n_tot": [60e6]}}'},
                                {"out_bin" : '{"radii": {"moms": [0], "drwt": "wet", "nbin": 26, "lnli": "log", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"aqq": 1, "rght": 0.0001, "moms": [0], "drwt": "wet", "nbin": 26, "lnli": "log", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": "aqq", "moms": [0], "drwt": "wet", "nbin": 26, "lnli": "log", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0001, "moms": [0], "drwt": "wet", "nbin": 26, "lnli": "log", "left": "aqq"}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0, "moms": [0], "drwt": "wet", "nbin": 26, "lnli": "log", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0001, "moms": [0], "drwt": "aqq", "nbin": 26, "lnli": "log", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0001, "moms": [0], "drwt": "wet", "nbin": 26, "lnli": "aqq", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0001, "moms": [0], "drwt": "wet", "nbin": "aqq", "lnli": "log", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0001, "moms": "aqq", "drwt": "wet", "nbin": 26, "lnli": "aqq", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0001, "moms": ["aqq"], "drwt": "wet", "nbin": 26, "lnli": "aqq", "left": 1e-09}}'},
                                {"out_bin" : '{"radii": {"rght": 0.0001, "moms": [0], "drwt": "wet", "nbin": 26, "lnli": "aqq", "left": 1e-09}}'}
                             ])

def test_args(tmpdir, arg):
    """checking if parcel rises exceptions when arguments don't make sense """
    str_f = str(tmpdir.join("test_pcl.nc"))
    with pytest.raises(Exception) as excinfo:
        pc.parcel(outfile=str_f,  **arg)
