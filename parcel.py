#!/usr/bin/env python

from argparse import ArgumentParser, RawTextHelpFormatter

from libcloudphxx import common, lgrngn
from libcloudphxx import git_revision as libcloud_version

from distutils.version import StrictVersion
from scipy import __version__ as scipy_version
assert StrictVersion(scipy_version) >= StrictVersion("0.13"), "see https://github.com/scipy/scipy/pull/491"

from scipy.io import netcdf
import inspect, numpy as np
import pdb

import subprocess
parcel_version = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip()

_Chem_ga_id = ["SO2", "H2O2", "O3"]
_Chem_aq_id = _Chem_ga_id + ["HSO3"]

_Chem_id = {
  "SO2"  : lgrngn.chem_species_t.SO2,
  "H2O2" : lgrngn.chem_species_t.H2O2,
  "O3"   : lgrngn.chem_species_t.O3,
  "HSO3" : lgrngn.chem_species_t.HSO3                                                      
}

def _micro_init(opts, state, info):
  # sanity check
  _stats(state, info)
  if (state["RH"] > 1): raise Exception("Please supply initial T,p,r_v below supersaturation")

  # using nested function to get access to opts
  def lognormal(lnr):
    from math import exp, log, sqrt, pi
    return opts["n_tot"] * exp(
      -(lnr - log(opts["mean_r"]))**2 / 2 / log(opts["gstdev"])**2
    ) / log(opts["gstdev"]) / sqrt(2*pi);

  # lagrangian scheme options
  opts_init = lgrngn.opts_init_t()  
  for opt in ["dt",]:  
    setattr(opts_init, opt, opts[opt])
  opts_init.sd_conc_mean = opts["sd_conc"]
  opts_init.dry_distros = {opts["kappa"]:lognormal}
  opts_init.kernel = lgrngn.kernel_t.geometric #TODO: will not be needed soon (libcloud PR #89)
  opts_init.chem_switch = True 

  # initialitation
  micro = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  micro.init(state["th_d"], state["r_v"], state["rhod"])
  return micro

def _micro_step(micro, state, info, chem_gas):
  libopts = lgrngn.opts_t()
  for id in chem_gas:
    libopts.chem_gas[_Chem_id[id]] = chem_gas[id]

  micro.step_sync(libopts, state["th_d"], state["r_v"], state["rhod"]) 

  # new = diag_chem(_Chem_id[id])
  for id in chem_gas:
    old = 0 #libopts.chem_gas[_Chem_id[id]]
    
    new = 0 #np.frombuffer(micro.outbuf())    
    chem_gas[id] -= (new - old)

def _stats(state, info):
  state["T"] = np.array([common.T(state["th_d"][0], state["rhod"][0])])
  state["RH"] = state["p"] * state["r_v"] / (state["r_v"] + common.eps) / common.p_vs(state["T"][0])
  info["RH_max"] = max(info["RH_max"], state["RH"])

def _histo(bins, micro, opts, chem_aq):
  r_min = 0
  i = 0
  for r_max in opts["radii"]:
    micro.diag_wet_rng(r_min, r_max)

    micro.diag_wet_mom(0) # #/kg dry air
    bins["conc"][i] = np.frombuffer(micro.outbuf())

    for id in _Chem_aq_id:
      micro.diag_chem(_Chem_id[id])
      chem_aq[id][i] = np.frombuffer(micro.outbuf())

    r_min = r_max
    i += 1

def _output_init(opts):
  # file & dimensions
  fout = netcdf.netcdf_file(opts["outfile"], 'w')
  fout.createDimension('t', None)
  fout.createDimension('radii', opts["radii"].shape[0]) #TODO: r_d, cloud only; #TODO: r_w vs. r_v - might be misleading
  
  units = {"z" : "m", "t" : "s", "r_v" : "kg/kg", "th_d" : "K", "rhod" : "kg/m3", 
    "p" : "Pa", "T" : "K", "RH" : "1", "conc" : "(kg of dry air)^-1"
  }
  for id in _Chem_id:
    units[id] = "todo"

  for name, unit in units.iteritems():
    if name in _Chem_aq_id + ["conc"]:
      dims = ('t','radii')
    else:
      dims = ('t',)

    fout.createVariable(name, 'd', dims)
    fout.variables[name].unit = unit

  return fout

def _output_save(fout, state, rec):
  for var, val in state.iteritems():
    fout.variables[var][rec] = val

def _save_attrs(fout, dictnr):
  for var, val in dictnr.iteritems():
    setattr(fout, var, val)

def _output(fout, opts, micro, bins, state, chem_gas, chem_aq, rec):
  _histo(bins, micro, opts, chem_aq)
  _output_save(fout, state, rec)
  _output_save(fout, bins, rec)
  _output_save(fout, chem_aq, rec) 
  _output_save(fout, chem_gas, rec)

 
def parcel(dt=.1, z_max=200., w=1., T_0=300., p_0=101300., r_0=.022, outfile="test.nc", 
  outfreq=100, sd_conc=64., kappa=.5,
  mean_r = .04e-6 / 2, gstdev  = 1.4, n_tot  = 60.e6, 
  radii = 1e-6 * pow(10, -3 + np.arange(26) * .2), 
  SO2_0 = 44., O3_0 = 44., H2O2_0 = 44.
):
  """
  Args:
    dt      (Optional[float]):   timestep [s]
    z_max   (Optional[float]):   maximum vertical displacement [m]
    w       (Optional[float]):   updraft velocity [m/s]
    T_0     (Optional[float]):   initial temperature [K]
    p_0     (Optional[float]):   initial pressure [Pa]
    r_0     (Optional[float]):   initial water vapour mass mixing ratio [kg/kg]
    outfile (Optional[string]):  output netCDF file name
    outfreq (Optional[int]):     output interval (in number of time steps)
    sd_conc (Optional[int]):     number of moving bins (super-droplets)
    kappa   (Optional[float]):   kappa hygroscopicity parameter (see doi:10.5194/acp-7-1961-2007)
    mean_r  (Optional[float]):   lognormal distribution mode diameter [m]
    gstdev  (Optional[float]):   lognormal distribution geometric standard deviation [1]
    n_tot   (Optional[float]):   lognormal distribution total concentration under standard 
                                 conditions (T=20C, p=1013.25 hPa, rv=0) [m-3]
    radii   (Optional[ndarray]): right bin edges for spectrum output [m]
                                 (left edge of the first bin equals 0)
    SO2_0   (Optional[float]):   initial SO2 TODO [TODO]
    O3_0    (Optional[float]):   initial O3 TODO [TODO]
    H2O2_0  (Optional[float]):   initial H2O2 TODO [TODO]
  """
  _arguments_checking(locals())

  # packing function arguments into "opts" dictionary
  args, _, _, _ = inspect.getargvalues(inspect.currentframe())
  opts = dict(zip(args, [locals()[k] for k in args]))

  th_0 = T_0 * (common.p_1000 / p_0)**(common.R_d / common.c_pd)
  nt = int(z_max / (w * dt))
  state = {
    "t" : 0, "z" : 0,
    "r_v" : np.array([r_0]), "p" : p_0,
    "th_d" : np.array([common.th_std2dry(th_0, r_0)]), 
    "rhod" : np.array([common.rhod(p_0, th_0, r_0)]),
    "T" : None, "RH" : None
  }
  info = { "RH_max" : 0, "libcloud_Git_revision" : libcloud_version, 
           "parcel_Git_revision" : parcel_version }
  bins = { "conc" : np.empty((radii.shape[0],)) }
  chem_gas = { "SO2" : SO2_0, "O3" : O3_0, "H2O2" : H2O2_0 }
  chem_aq = dict(zip(_Chem_aq_id, len(_Chem_aq_id)*[np.empty(radii.shape[0])]))
  with _output_init(opts) as fout:
    # t=0 : init & save
    micro = _micro_init(opts, state, info)
    _output(fout, opts, micro, bins, state, chem_gas, chem_aq, 0)

    # timestepping
    for it in range(1,nt+1):
      # diagnostics
      # the reasons to use analytic solution:
      # - independent of dt
      # - same as in 2D kinematic model
      state["z"] += w * dt
      state["t"] = it * dt
      state["p"] = common.p_hydro(state["z"], th_0, r_0, 0, p_0)
      state["rhod"][0] = common.rhod(state["p"], th_0, r_0)

      # microphysics
      _micro_step(micro, state, info, chem_gas)
      _stats(state, info)
    
      # TODO: only if user wants to stop @ RH_max
      #if (state["RH"] < info["RH_max"]): break

      # output
      if (it % outfreq == 0): 
        rec = it/outfreq
        _output(fout, opts, micro, bins, state, chem_gas, chem_aq, rec)
 
    _save_attrs(fout, info)
    _save_attrs(fout, opts)

    
def _arguments_checking(args):
  if (args["gstdev"] == 1): raise Exception("standar deviation should be != 1 to avoid monodisperse distribution")
  if (args["T_0"] < 273.15): raise Exception("temperature should be larger than 0C - microphysics works only for warm clouds")
  if (args["r_0"] < 0): raise Exception("water vapour should be larger than 0")
  if (args["w"] < 0): raise Exception("vertical velocity should be larger than 0")



# ensuring that pure "import parcel" does not trigger any simulation
if __name__ == '__main__':

  # getting list of argument names and their default values
  name, _, _, dflt = inspect.getargspec(parcel)
  opts = dict(zip(name[-len(dflt):], dflt))

  # handling all parcel() arguments as command-line arguments
  prsr = ArgumentParser(add_help=True, description=parcel.__doc__, formatter_class=RawTextHelpFormatter)
  for k in opts:
    prsr.add_argument('--' + k, 
      default=opts[k], 
      help = "(default: %(default)s)",
      # reading in ndarrays as lists (see comment below ****)
      type = (type(opts[k]) if type(opts[k]) != np.ndarray else type(opts[k][0])),
      nargs = ('?'          if type(opts[k]) != np.ndarray else '+')
    )
  args = vars(prsr.parse_args())

  # converting lists into ndarrays (see comment abowe ****)
  for k in opts:
    if type(opts[k]) == np.ndarray:
      args[k] = np.fromiter(args[k], dtype=opts[k].dtype)

  # executing parcel() with command-line arguments unpacked - treated as keyword arguments 
  parcel(**args)
