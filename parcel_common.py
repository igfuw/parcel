#!/usr/bin/env python
import numpy as np
from libcloudphxx import common, lgrngn

# id_str     id_int (gas phase chemistry labels)
_Chem_g_id = {
  "SO2_g"  : lgrngn.chem_species_t.SO2,
  "H2O2_g" : lgrngn.chem_species_t.H2O2,
  "O3_g"   : lgrngn.chem_species_t.O3,
  "HNO3_g" : lgrngn.chem_species_t.HNO3,
  "NH3_g"  : lgrngn.chem_species_t.NH3,
  "CO2_g"  : lgrngn.chem_species_t.CO2
}

# id_str     id_int (aqueous phase chemistry labels)
_Chem_a_id = {
  "SO2_a"  : lgrngn.chem_species_t.SO2,
  "H2O2_a" : lgrngn.chem_species_t.H2O2,
  "O3_a"   : lgrngn.chem_species_t.O3,
  "CO2_a"  : lgrngn.chem_species_t.CO2,
  "HNO3_a" : lgrngn.chem_species_t.HNO3,
  "NH3_a"  : lgrngn.chem_species_t.NH3,
  "H"      : lgrngn.chem_species_t.H,
  "S_VI"   : lgrngn.chem_species_t.S_VI
}

def _w_eval(w, t):
  """Evaluate vertical velocity at time t.

  Supported forms for `w`:
    - float/int: constant vertical velocity
    - callable:  w(t) -> float
    - str:       Python expression in variable `t` (seconds), e.g. "1 + 0.5*np.sin(2*np.pi*t/60)"
                Available names: t, np

  Note: expression strings are `eval`'d with a restricted global namespace.
  """
  if callable(w):
    return float(w(t))
  if isinstance(w, (int, float, np.floating)):
    return float(w)
  if isinstance(w, str):
    # restricted eval environment
    return float(eval(w, {"__builtins__": {}}, {"t": float(t), "np": np}))
  raise TypeError("w must be a number, a callable w(t), or an expression string")


class lognormal(object):
  def __init__(self, mean_r, gstdev, n_tot):
    self.mean_r = mean_r
    self.gstdev = gstdev
    self.n_tot = n_tot

  def __call__(self, lnr):
    from math import exp, log, sqrt, pi
    return self.n_tot * exp(
      -(lnr - log(self.mean_r))**2 / 2 / log(self.gstdev)**2
    ) / log(self.gstdev) / sqrt(2*pi);

class sum_of_lognormals(object):
  def __init__(self, lognormals=[]):
    self.lognormals = lognormals

  def __call__(self, lnr):
    res = 0.
    for lognormal in self.lognormals:
      res += lognormal(lnr)
    return res


def _stats(state, info):
  state["T"] = np.array([common.T(state["th_d"][0], state["rhod"][0])])
  state["RH"] = state["p"] * state["r_v"] / (state["r_v"] + common.eps) / common.p_vs(state["T"][0])
  state["T_blk"] = np.array([state["th_d"][0] * common.exner(state["p"])])
  state["RH_blk"] = state["r_v"] / common.r_vs(state["T_blk"][0], state["p"])
  info["RH_max"] = max(info["RH_max"], state["RH"])


def _p_hydro_const_rho(dz, p, rho):
  # hydrostatic pressure assuming constatnt density
  return p - rho * common.g * dz


def _p_hydro_const_th_rv(z_lev, p_0, th_std, r_v, z_0=0.):
  # hydrostatic pressure assuming constatnt theta and r_v
  return common.p_hydro(z_lev, th_std, r_v, z_0, p_0)

def _arguments_checking(opts, spectra, aerosol, ice_switch):
  if opts["T_0"] < 273.15 and ice_switch == False:
    raise Exception("temperature should be larger than 0C if ice_switch=False")
  elif ((opts["r_0"] >= 0) and (opts["RH_0"] >= 0)):
    raise Exception("both r_0 and RH_0 specified, please use only one")
  # if opts["w"] < 0:
  #   raise Exception("vertical velocity should be larger than 0")
  if opts["t"] is not None and opts["z_max"] is not None:
    raise ValueError("Specify only one stop condition: either z_max or t (not both)")
  if opts["t"] is None and opts["z_max"] is None:
    raise ValueError("You must specify exactly one stop condition: z_max or t")
  if opts["t"] is not None and opts["t"] <= 0:
    raise ValueError("t must be > 0")
  if opts["z_max"] is not None and opts["z_max"] <= 0:
    raise ValueError("z_max must be > 0")
  w0 = _w_eval(opts["w"], 0.0)
  if w0 <= 0 and isinstance(opts["w"], (int, float, np.floating)) and opts["z_max"] is not None:
    raise ValueError("For constant w with z_max stop, expected w>0 to reach z_max")

  if aerosol is not None:
    for name, dct in aerosol.items():
      # TODO: check if name is valid netCDF identifier
      # (http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM/Identifiers.html)
      keys = ["kappa", "mean_r", "n_tot", "gstdev"]
      for key in keys:
        if key not in dct:
          raise Exception(">>" + key + "<< is missing in aerosol[" + name + "]")
      for key in dct:
        if key not in keys:
          raise Exception("invalid key >>" + key + "<< in aerosol[" + name + "]")
      if dct["kappa"] <= 0:
        raise Exception("kappa hygroscopicity parameter should be larger than 0 for aerosol[" + name + "]")
      if type(dct["mean_r"]) != list:
          raise Exception(">>mean_r<< key in aerosol["+ name +"] must be a list")
      if type(dct["gstdev"]) != list:
          raise Exception(">>gstdev<< key in aerosol["+ name +"] must be a list")
      if type(dct["n_tot"]) != list:
          raise Exception(">>n_tot<< key in aerosol["+ name +"] must be a list")
      if not len(dct["mean_r"]) == len(dct["n_tot"]) == len(dct["gstdev"]):
        raise Exception("mean_r, n_tot and gstdev lists should have same sizes for aerosol[" + name + "]")
      for mean_r in dct["mean_r"]:
        if mean_r <= 0:
          raise Exception("mean radius should be > 0 for aerosol[" + name + "]")
      for n_tot in dct["n_tot"]:
        if n_tot <= 0:
          raise Exception("concentration should be > 0 for aerosol[" + name + "]")
      for gstdev in dct["gstdev"]:
        if gstdev <= 0:
          raise Exception("standard deviation should be > 0 for aerosol[" + name + "]")
      # necessary?
        if gstdev == 1.:
          raise Exception("standard deviation should be != 1 to avoid monodisperse distribution for aerosol[" + name + "]")

  for name, dct in spectra.items():
    # TODO: check if name is valid netCDF identifier
    # (http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM/Identifiers.html)
    keys = ["left", "rght", "nbin", "drwt", "lnli", "moms"]
    for key in keys:
      if key not in dct:
        raise Exception(">>" + key + "<< is missing in out_bin[" + name + "]")
    for key in dct:
      if key not in keys:
        raise Exception("invalid key >>" + key + "<< in out_bin[" + name + "]")
    if type(dct["left"]) not in [int, float]:
        raise Exception(">>left<< in out_bin["+ name +"] must be int or float")
    if type(dct["rght"]) not in [int, float]:
        raise Exception(">>rght<< in out_bin["+ name +"] must be int or float")
    if dct["left"] >= dct["rght"]:
        raise Exception(">>left<< is greater than >>rght<< in out_bin["+ name +"]")
    if dct["drwt"] not in ["dry", "wet"]:
        raise Exception(">>drwt<< key in out_bin["+ name +"] must be either >>dry<< or >>wet<<")
    if dct["lnli"] not in ["lin", "log"]:
        raise Exception(">>lnli<< key in out_bin["+ name +"] must be either >>lin<< or >>log<<")
    if type(dct["nbin"]) != int:
        raise Exception(">>nbin<< key in out_bin["+ name +"] must be an integer number")
    if type(dct["moms"]) != list:
        raise Exception(">>moms<< key in out_bin["+ name +"] must be a list")
    for mom in dct["moms"]:
        if (type(mom) != int):
          if (mom not in list(_Chem_a_id.keys())):
            raise Exception(">>moms<< key in out_bin["+ name +"] must be a list of integer numbers or valid chemical compounds (" +str(list(_Chem_a_id.keys())) + ")")

def _init_sanity_check(state, info):
    _stats(state, info)
    if (state["RH"] > 1): 
        raise Exception("Please supply initial T,p,r_v below supersaturation")
