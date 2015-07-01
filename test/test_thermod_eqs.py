import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
import parcel as pc
from scipy.io import netcdf
import pytest


@pytest.mark.parametrize("rho_cons, th_d, r_v", [(1., 300., 0), (1.3, 300., 1.e-4), 
                                                 (1.1, 280., 1.e-5), (.8, 300., 1.e-5)])
def test_press_compat(rho_cons, th_d, r_v, dz=20, p0=8.e4):
    """ checking if thermodynamic eqs. used to calculate pressure and density are consistent""" 
    p1 = pc._p_hydro_const_rho(dz, p0, rho_cons)
    rhod = pc.common.rhod(p1, pc.common.th_dry2std(th_d, r_v), r_v)
    p2 = pc.common.p(rhod, r_v, pc.common.T(th_d, rhod))
    assert abs(p1 - p2) < 1.e-9 * p1 
