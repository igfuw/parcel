#!/usr/bin/env python

from argparse import ArgumentParser, RawTextHelpFormatter

from libcloudphxx import common, lgrngn
from libcloudphxx import git_revision as libcloud_version

from distutils.version import StrictVersion
from scipy import __version__ as scipy_version
assert StrictVersion(scipy_version) >= StrictVersion("0.13"), "see https://github.com/scipy/scipy/pull/491"

from scipy.io import netcdf
import re, inspect, numpy as np
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

def _output_bins(fout, t, micro, opts):
  for dim, nbin in fout.dimensions.iteritems():
    if (dim == 't'): continue
    for b in range(nbin):
      micro.diag_wet_rng(
	fout.variables[dim+"_rl"][b],
	fout.variables[dim+"_rl"][b] + fout.variables[dim+"_dr"][b]
      )
      for v in fout.variables.iterkeys():
        if v.startswith(dim+"_"): 
          #TODO: chemia
          rgxp = re.search('^'+dim+'_m(\d+)$', v)
          if rgxp != None:
            m = int(rgxp.groups()[0])
            micro.diag_wet_mom(m)
            fout.variables[dim+'_m'+str(m)][t, b] = np.frombuffer(micro.outbuf())

def _output_init(opts):
  # file & dimensions
  fout = netcdf.netcdf_file(opts["outfile"], 'w')
  fout.createDimension('t', None)
  for e in opts["out_wet"]:
    (
      name   ,left      ,rght      ,nbin   ,lnli ,moms 
    ) = [t(s) for t,s in zip((
      str    ,float     ,float     ,int    ,str  ,str
    ),re.search(
      '^(\w+):([\d.e-]+)/([\d.e-]+)/([\d]+)/(\w+)/([\d,]+)$', 
      e
    ).groups())]
    fout.createDimension(name, nbin) 
    fout.createVariable(name+'_rl', 'd', (name,))
    fout.createVariable(name+'_dr', 'd', (name,))
    if lnli == 'log':
      from math import exp, log
      dlnr = (log(rght) - log(left)) / nbin
      allbins = np.exp(log(left) + np.arange(nbin+1) * dlnr)
      fout.variables[name+'_rl'][:] = allbins[0:-1]
      fout.variables[name+'_dr'][:] = allbins[1:] - allbins[0:-1]
    elif lnli == 'lin':
      dr = (rght - left) / nbin
      fout.variables[name+'_rl'][:] = left * np.arange(nbin) * dr
      fout.variables[name+'_dr'][:] = dr
    else:
      raise exception('scale type can be either log or lin')
    for m in moms.split(','):
      #TODO
      #if (m in _Chem_aq_id):
      #	fout.createVariable(name+'_'+m, 'd', ('t',name))
      #	fout.variables[name+'_m'+m].unit = 'kg of chem species dissolved in cloud droplets (kg of dry air)^-1'
      #else
        assert(str(int(m))==m)
	fout.createVariable(name+'_m'+m, 'd', ('t',name))
	fout.variables[name+'_m'+m].unit = 'm^'+m+' (kg of dry air)^-1'
  
  units = {"z" : "m", "t" : "s", "r_v" : "kg/kg", "th_d" : "K", "rhod" : "kg/m3", 
    "p" : "Pa", "T" : "K", "RH" : "1"
  }

  # TODO: 
  # for id_str in _Chem_g_id.iterkeys():
  #   units[id_str] = "gas volume concentration (mole fraction) [1]"

  for name, unit in units.iteritems():
    fout.createVariable(name, 'd', ('t',))
    fout.variables[name].unit = unit

  return fout

def _output_save(fout, state, rec):
  for var, val in state.iteritems():
    fout.variables[var][rec] = val

def _save_attrs(fout, dictnr):
  for var, val in dictnr.iteritems():
    setattr(fout, var, val)

def _output(fout, opts, micro, state, chem_gas, rec):
  _output_bins(fout, rec, micro, opts)
  _output_save(fout, state, rec)
  _output_save(fout, chem_gas, rec)

def _p_hydro_const_rho(dz, p, rho):
  # hydrostatic pressure assuming constatnt density
  return p - rho * common.g * dz

def _p_hydro_const_th_rv(z_lev, p_0, th_std, r_v, z_0=0.):
  # hydrostatic pressure assuming constatnt theta and r_v
  return common.p_hydro(z_lev, th_std, r_v, z_0, p_0)
 
def parcel(dt=.1, z_max=200., w=1., T_0=300., p_0=101300., r_0=.022, 
  outfile="test.nc", 
  pprof="pprof_piecewise_const_rhod",
  outfreq=100, sd_conc=64., kappa=.5,
  mean_r = .04e-6 / 2, gstdev  = 1.4, n_tot  = 60.e6, 
  out_wet = ["radii:1e-9/1e-4/26/log/0"], 
  #radii = 1e-6 * pow(10, -3 + np.arange(26) * .2), 
  SO2_0 = 0., O3_0 = 0., H2O2_0 = 0.
):
  """
  Args:
    dt      (Optional[float]):    timestep [s]
    z_max   (Optional[float]):    maximum vertical displacement [m]
    w       (Optional[float]):    updraft velocity [m/s]
    T_0     (Optional[float]):    initial temperature [K]
    p_0     (Optional[float]):    initial pressure [Pa]
    r_0     (Optional[float]):    initial water vapour mass mixing ratio [kg/kg]
    outfile (Optional[string]):   output netCDF file name
    outfreq (Optional[int]):      output interval (in number of time steps)
    sd_conc (Optional[int]):      number of moving bins (super-droplets)
    kappa   (Optional[float]):    kappa hygroscopicity parameter (see doi:10.5194/acp-7-1961-2007)
    mean_r  (Optional[float]):    lognormal distribution mode diameter [m]
    gstdev  (Optional[float]):    lognormal distribution geometric standard deviation [1]
    n_tot   (Optional[float]):    lognormal distribution total concentration under standard 
                                  conditions (T=20C, p=1013.25 hPa, rv=0) [m-3]
    out_wet (Optional[str list]): array of strings defining spectrum diagnostics, e.g.:
                                  ["radii:0/1e-4/26/log/0","cloud:.5e-6/25e-6/49/lin/0,1,2,3"]
                                  will generate five output spectra:
                                  - 0-th spectrum moment for 26 bins spaced logarithmically between 0 and 1e-4 m
                                  - 0,1,2 & 3-rd moments for 49 bins spaced linearly between .5e-6 and 25e-6
    SO2_0   (Optional[float]):    initial SO2 TODO [TODO]
    O3_0    (Optional[float]):    initial O3 TODO [TODO]
    H2O2_0  (Optional[float]):    initial H2O2 TODO [TODO]
    pprof   (Optional[string]):   method to calculate pressure profile used to calculate 
                                  dry air density that is used by the super-droplet scheme
                                  valid options are: pprof_const_th_rv, pprof_const_rhod, pprof_piecewise_const_rhod
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
  #bins = { "conc" : np.empty((radii.shape[0],)) }
  chem_gas = { }#"SO2" : SO2_0, "O3" : O3_0, "H2O2" : H2O2_0 }
  with _output_init(opts) as fout:
    # t=0 : init & save
    micro = _micro_init(opts, state, info)
    _output(fout, opts, micro, state, chem_gas, 0)

    # timestepping
    for it in range(1,nt+1):
      # diagnostics
      # the reasons to use analytic solution:
      # - independent of dt
      # - same as in 2D kinematic model
      state["z"] += w * dt
      state["t"] = it * dt

      # pressure
      if pprof == "pprof_const_th_rv":
        # as in icicle model
        p_hydro = _p_hydro_const_th_rv(state["z"], p_0, th_0, r_0)
      elif pprof == "pprof_const_rhod":
        # as in Grabowski and Wang 2009
        rho = 1.13 # kg/m3  1.13 
        state["p"] = _p_hydro_const_rho(state["z"], p_0, rho) 

      elif pprof == "pprof_piecewise_const_rhod":
        # as in Grabowski and Wang 2009 but calculating pressure
        # for rho piecewise constant per each time step
        state["p"] = _p_hydro_const_rho(w*dt, state["p"], state["rhod"][0])

      else: assert(False)

      # dry air density
      if pprof == "pprof_const_th_rv":
        state["rhod"][0] = common.rhod(p_hydro, th_0, r_0)
        state["p"] = common.p(
          state["rhod"][0],
          state["r_v"][0],
          common.T(state["th_d"][0], state["rhod"][0])
        )

      else:
        state["rhod"][0] = common.rhod(
          state["p"], 
          common.th_dry2std(state["th_d"][0], state["r_v"][0]), 
          state["r_v"][0]
        )

      # microphysics
      _micro_step(micro, state, info, chem_gas)
      _stats(state, info)
    
      # TODO: only if user wants to stop @ RH_max
      #if (state["RH"] < info["RH_max"]): break

      # output
      if (it % outfreq == 0): 
        rec = it/outfreq
        _output(fout, opts, micro, state, chem_gas, rec)
 
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
      type = (type(opts[k]) if type(opts[k]) != list else type(opts[k][0])),
      nargs = ('?'          if type(opts[k]) != list else '+')
    )
  args = vars(prsr.parse_args())

  # executing parcel() with command-line arguments unpacked - treated as keyword arguments 
  parcel(**args)
