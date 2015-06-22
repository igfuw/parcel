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

# id_str     id_int
_Chem_g_id = {
  "SO2_g"  : lgrngn.chem_species_t.SO2, 
  "H2O2_g" : lgrngn.chem_species_t.H2O2, 
  "O3_g"   : lgrngn.chem_species_t.O3
}

# id_str     id_int
_Chem_a_id = {
  "SO2_a"  : lgrngn.chem_species_t.SO2, 
  "H2O2_a" : lgrngn.chem_species_t.H2O2, 
  "O3_a"   : lgrngn.chem_species_t.O3,
  "HSO3_a" : lgrngn.chem_species_t.HSO3 
}

# id_int   ...
_molar_mass = { #... TODO molar mass of dissoved one -> + M_H2O
  lgrngn.chem_species_t.SO2  : common.M_SO2,
  lgrngn.chem_species_t.H2O2 : common.M_H2O2,
  lgrngn.chem_species_t.O3   : common.M_O3
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
  opts_init.chem_rho = opts["chem_rho"]

  # switching off chemistry if all initial volume conc. equal zero
  opts_init.chem_switch = False
  for id_str in _Chem_g_id.iterkeys():
    if opts[id_str + "_0"] != 0: opts_init.chem_switch = True

  # initialisation
  micro = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  micro.init(state["th_d"], state["r_v"], state["rhod"])
  return micro

def _micro_step(micro, state, info, opts):
  libopts = lgrngn.opts_t()
  libopts.cond = True
  libopts.chem = micro.opts_init.chem_switch
  libopts.coal = False
  libopts.adve = False
  libopts.sedi = False

  if micro.opts_init.chem_switch:
    tmp = {}
    for id_str, id_int in _Chem_g_id.iteritems():
      tmp[id_int] = state[id_str]
    libopts.chem_gas = tmp

  print "old rv = ", state["r_v"]
  micro.step_sync(libopts, state["th_d"], state["r_v"], state["rhod"]) 
  print "new rv = ", state["r_v"]
  print " "

  micro.step_async(libopts)
  _stats(state, info) # give updated T needed for chemistry below

  if micro.opts_init.chem_switch:
    micro.diag_all() # selecting all particles
    for id_str, id_int in _Chem_g_id.iteritems():
      if opts['chem_sys'] == 'closed':

        old = state[id_str.replace('_g', '_a')]

        micro.diag_chem(id_int)
        new = np.frombuffer(micro.outbuf()) 
      
        # since p & rhod are the "new" ones, for consistency we also use new T (_stats called above)
        state[id_str] -= (new[0] - old) * state["rhod"][0] * common.R * state["T"][0] / _molar_mass[id_int] / state["p"]
        state[id_str.replace('_g', '_a')] = new[0]
        #print "state[ ", id_str.replace('_g', '_a'), " ] = ", np.frombuffer(micro.outbuf())[0]
      elif opts['chem_sys'] == 'open':
        micro.diag_chem(id_int)
        state[id_str.replace('_g', '_a')] = np.frombuffer(micro.outbuf())[0]
      else:
        raise exception(
          "Expected chem_sys options are: 'open', 'closed'."
          "Type: help(parcel) for more help."
        )
 
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
    bins["conc_bins"][i] = np.frombuffer(micro.outbuf())

    if micro.opts_init.chem_switch:
      for id_str, id_int in _Chem_a_id.iteritems():
        micro.diag_chem(id_int)
        chem_aq[id_str + '_bins'][i] = np.frombuffer(micro.outbuf())
    
    r_min = r_max
    i += 1

def _output_init(micro, opts):
  # file & dimensions
  fout = netcdf.netcdf_file(opts["outfile"], 'w')
  fout.createDimension('t', None)
  fout.createDimension('radii', opts["radii"].shape[0]) #TODO: r_d, cloud only; #TODO: r_w vs. r_v - might be misleading
  
  units = {"z" : "m", "t" : "s", "r_v" : "kg/kg", "th_d" : "K", "rhod" : "kg/m3", 
    "p" : "Pa", "T" : "K", "RH" : "1", "conc_bins" : "(kg of dry air)^-1" 
  }
  if micro.opts_init.chem_switch:
    for id_str in _Chem_a_id.iterkeys():
      unit = "kilograms of chem species dissolved within droplets / kilograms of dry air"
      units[id_str] = unit 
      units[id_str +  '_bins'] = unit 
    for id_str in _Chem_g_id.iterkeys():
      units[id_str] = "gas volume concentration (mole fraction) [1]"

  for name, unit in units.iteritems():
    if name.endswith('_bins'):
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

def _output(fout, opts, micro, bins, state, chem_aq, rec):
  _histo(bins, micro, opts, chem_aq)
  _output_save(fout, state, rec)
  _output_save(fout, bins, rec)
  if micro.opts_init.chem_switch:
    _output_save(fout, chem_aq, rec) 
 
def parcel(dt=1., z_max=200., w=1., T_0=300., p_0=101300., r_0=.022, outfile="test.nc", 
  outfreq=10, sd_conc=64., kappa=.5,
  mean_r = .04e-6 / 2, gstdev  = 1.4, n_tot  = 60.e6, 
  radii = 1e-6 * pow(10, -3 + np.arange(26) * .2), 
#  SO2_g_0 = 0., O3_g_0 = 0., H2O2_g_0 = 0.,
  SO2_g_0 = 200e-12, O3_g_0 = 50e-9, H2O2_g_0 = 500e-12,
  chem_sys = 'closed', 
  chem_rho = 1.8e-3
):
  """
  Args:
    dt       (Optional[float]):   timestep [s]
    z_max    (Optional[float]):   maximum vertical displacement [m]
    w        (Optional[float]):   updraft velocity [m/s]
    T_0      (Optional[float]):   initial temperature [K]
    p_0      (Optional[float]):   initial pressure [Pa]
    r_0      (Optional[float]):   initial water vapour mass mixing ratio [kg/kg]
    outfile  (Optional[string]):  output netCDF file name
    outfreq  (Optional[int]):     output interval (in number of time steps)
    sd_conc  (Optional[int]):     number of moving bins (super-droplets)
    kappa    (Optional[float]):   kappa hygroscopicity parameter (see doi:10.5194/acp-7-1961-2007)
    mean_r   (Optional[float]):   lognormal distribution mode diameter [m]
    gstdev   (Optional[float]):   lognormal distribution geometric standard deviation [1]
    n_tot    (Optional[float]):   lognormal distribution total concentration under standard 
                                  conditions (T=20C, p=1013.25 hPa, rv=0) [m-3]
    radii    (Optional[ndarray]): right bin edges for spectrum output [m]
                                  (left edge of the first bin equals 0)
    SO2_g_0  (Optional[float]):   initial SO2  gas volume concentration (mole fraction) [1]
    O3_g_0   (Optional[float]):   initial O3   gas volume concentration (mole fraction) [1]
    H2O2_g_0 (Optional[float]):   initial H2O2 gas volume concentration (mole fraction) [1]
    chem_sys (Optional[string]):  accepted values: 'open' or 'closed'
                                  (in open/closed system gas volume concentration in the air doesn't/does change 
                                   due to chemical reactions)
  """
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
  bins = { "conc_bins" : np.empty((radii.shape[0],)) }
  tmp = _Chem_a_id.keys()
  chem_aq = dict(zip([x + '_bins' for x in tmp], len(tmp)*[np.empty(radii.shape[0])]))

  micro = _micro_init(opts, state, info)
  with _output_init(micro, opts) as fout:
    # adding chem state vars
    if micro.opts_init.chem_switch:
      state.update({ "SO2_g" : SO2_g_0, "O3_g" : O3_g_0, "H2O2_g" : H2O2_g_0 })
      state.update({ "SO2_a" : 0.,      "O3_a" : 0.,     "H2O2_a" : 0.      , "HSO3_a" : 0})

    # save state @ t=0
    _output(fout, opts, micro, bins, state, chem_aq, 0)

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
      _micro_step(micro, state, info, opts)
   
      # TODO: only if user wants to stop @ RH_max
      #if (state["RH"] < info["RH_max"]): break
 
      # output
      if (it % outfreq == 0): 
        rec = it/outfreq
        _output(fout, opts, micro, bins, state, chem_aq, rec)
 
    _save_attrs(fout, info)
    _save_attrs(fout, opts)


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

  # converting lists into ndarrays (see comment above ****)
  for k in opts:
    if type(opts[k]) == np.ndarray:
      args[k] = np.fromiter(args[k], dtype=opts[k].dtype)

  # executing parcel() with command-line arguments unpacked - treated as keyword arguments 
  parcel(**args)
