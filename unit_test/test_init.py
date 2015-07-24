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


@pytest.mark.parametrize("arg",[{"gstdev" : 1},
                                {"T_0" : 255},
                                {"r_0":-0.1},
                                {"w" : -1},
                                {"kappa" : 0.},
                                {"chem_sys" : "aqq"},
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
