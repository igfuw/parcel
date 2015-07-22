#!/usr/bin/env python

#import sys
#sys.path.insert(0, "../libcloudphxx/build/bindings/python/")

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
  "HSO3_a" : lgrngn.chem_species_t.HSO3,
  "SO3_a"  : lgrngn.chem_species_t.SO3,
  "HSO4_a" : lgrngn.chem_species_t.HSO4,
  "SO4_a"  : lgrngn.chem_species_t.SO4,
  "H"      : lgrngn.chem_species_t.H,
  "OH"     : lgrngn.chem_species_t.OH,
  "S_VI"   : lgrngn.chem_species_t.S_VI
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

  # switching on chemistry if either dissolving, dissociation or reactions are chosen
  opts_init.chem_switch = False
  if opts["chem_dsl"] or opts["chem_dsc"] or opts["chem_rct"]: opts_init.chem_switch = True

  # initialisation
  micro = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  micro.init(state["th_d"], state["r_v"], state["rhod"])
  return micro

def _micro_step(micro, state, info, opts, it):
  libopts = lgrngn.opts_t()
  libopts.cond = True
  libopts.coal = False
  libopts.adve = False
  libopts.sedi = False

  if micro.opts_init.chem_switch:
    tmp = {}
    for id_str, id_int in _Chem_g_id.iteritems():
      tmp[id_int] = state[id_str]
    libopts.chem_gas = tmp
    libopts.chem_dsl = opts["chem_dsl"]
    libopts.chem_dsc = opts["chem_dsc"]
    if it < opts["chem_spn"]:
        libopts.chem_rct = False
    else:
        libopts.chem_rct = opts["chem_rct"]

  micro.step_sync(libopts, state["th_d"], state["r_v"], state["rhod"]) 

  micro.step_async(libopts)
  _stats(state, info) # give updated T needed for chemistry below

  if micro.opts_init.chem_switch:
    micro.diag_all() # selecting all particles
    for id_str, id_int in _Chem_g_id.iteritems():
      if opts['chem_sys'] == 'closed':

        #TODO - some smarter way to add up decrease of gases due to dissolving into droplets
        # right now it's big - small all the way

        old = state[id_str.replace('_g', '_a')]

        micro.diag_chem(id_int)
        new = np.frombuffer(micro.outbuf()) 
      
        # since p & rhod are the "new" ones, for consistency we also use new T (_stats called above)
        state[id_str] -= (new[0] - old) * state["rhod"][0] * common.R * state["T"][0] / _molar_mass[id_int] / state["p"]
        assert state[id_str] >= 0

        state[id_str.replace('_g', '_a')] = new[0]

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

def _output_bins(fout, t, micro, opts):
  for dim, nbin in fout.dimensions.iteritems():
    if (dim == 't'): continue
    for b in range(nbin):
      drwt = 'wet' if dim+"_r_wet" in fout.variables else 'dry'
      if drwt == 'wet':
	micro.diag_wet_rng(
	  fout.variables[dim+"_r_wet"][b],
	  fout.variables[dim+"_r_wet"][b] + fout.variables[dim+"_dr_wet"][b]
	)
      elif drwt == 'dry':
	micro.diag_dry_rng(
	  fout.variables[dim+"_r_dry"][b],
	  fout.variables[dim+"_r_dry"][b] + fout.variables[dim+"_dr_dry"][b]
	)
      else: assert False
       
      for v in fout.variables.iterkeys():
        if v.startswith(dim+"_"):
          match = re.search('^'+dim+'_(\w+)$', v).groups()[0]
          if match in ['dr_wet', 'r_wet', 'dr_dry', 'r_dry']:
            pass
          elif match.startswith('m'):
            # calculating moments (they all have to start with m)
            mom = int(match[1:])
            if drwt == 'wet':
	      micro.diag_wet_mom(mom)
            elif drwt == 'dry':
	      micro.diag_dry_mom(mom)
            else: assert False
            fout.variables[dim+'_'+match][t, b] = np.frombuffer(micro.outbuf())
          else:
            # calculate chemistry
            micro.diag_chem(_Chem_a_id[match])
            fout.variables[dim+'_'+match][t, b] = np.frombuffer(micro.outbuf())
                      

def _output_init(micro, opts):
  # file & dimensions
  fout = netcdf.netcdf_file(opts["outfile"], 'w')
  fout.createDimension('t', None)
  for e in opts["out_bin"]:
    (
      name   ,left      ,rght      ,nbin   ,lnli ,drwt  ,moms 
    ) = [t(s) for t,s in zip((
      str    ,float     ,float     ,int    ,str  ,str   ,str
    ),re.search(
      '^(\w+):([\d.e-]+)/([\d.e-]+)/([\d]+)/(\w+)/(\w+)/([\d,\w]+)$', 
      e
    ).groups())]
    assert '_' not in name
    if drwt not in ['dry', 'wet']:
      raise exception('radius type can be either dry or wet')
    fout.createDimension(name, nbin) 

    tmp = name + '_r_' + drwt
    fout.createVariable(tmp, 'd', (name,))
    fout.variables[tmp].unit = "m"
    fout.variables[tmp].description = "particle wet radius (left bin edge)"

    tmp = name + '_dr_' + drwt
    fout.createVariable(tmp, 'd', (name,))
    fout.variables[tmp].unit = "m"
    fout.variables[tmp].description = "bin width"
    
    if lnli == 'log':
      from math import exp, log
      dlnr = (log(rght) - log(left)) / nbin
      allbins = np.exp(log(left) + np.arange(nbin+1) * dlnr)
      fout.variables[name+'_r_'+drwt][:] = allbins[0:-1]
      fout.variables[name+'_dr_'+drwt][:] = allbins[1:] - allbins[0:-1]
    elif lnli == 'lin':
      dr = (rght - left) / nbin
      fout.variables[name+'_r_'+drwt][:] = left + np.arange(nbin) * dr
      fout.variables[name+'_dr_'+drwt][:] = dr
    else:
      raise exception('scale type can be either log or lin')
    for m in moms.split(','):
      if (m in _Chem_a_id):
      	fout.createVariable(name+'_'+m, 'd', ('t',name))
      	fout.variables[name+'_'+m].unit = 'kg of chem species dissolved in cloud droplets (kg of dry air)^-1'
      else:
        assert(str(int(m))==m)
	fout.createVariable(name+'_m'+m, 'd', ('t',name))
	fout.variables[name+'_m'+m].unit = 'm^'+m+' (kg of dry air)^-1'
  
  units = {"z" : "m",  "t" : "s", "r_v" : "kg/kg", "th_d" : "K", "rhod" : "kg/m3", 
           "p" : "Pa", "T" : "K", "RH"  : "1"
  }

  if micro.opts_init.chem_switch:
    for id_str in _Chem_g_id.iterkeys():
      units[id_str] = "gas volume concentration (mole fraction) [1]"
      units[id_str.replace('_g', '_a')] = "kg of chem species dissolved in cloud droplets (kg of dry air)^-1"

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

def _output(fout, opts, micro, state, rec):
  _output_bins(fout, rec, micro, opts)
  _output_save(fout, state, rec)

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
  out_bin = ["radii:1e-9/1e-4/26/log/wet/0"], 
  SO2_g_0 = 0., O3_g_0 = 0., H2O2_g_0 = 0.,
  chem_sys = 'open',
  chem_dsl = False, chem_dsc = False, chem_rct = False, 
  chem_spn = 1,
  chem_rho = 1.8e3
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
    out_bin (Optional[str list]): array of strings defining spectrum diagnostics, e.g.:
                                  ["radii:0/1e-4/26/log/dry/0","cloud:.5e-6/25e-6/49/lin/wet/0,1,2,3"]
                                  will generate five output spectra:
                                  - 0-th spectrum moment for 26 bins spaced logarithmically between 0 and 1e-4 m 
                                    for dry radius
                                  - 0,1,2 & 3-rd moments for 49 bins spaced linearly between .5e-6 and 25e-6
                                    for wet radius
                                    (TODO - add chemistry output description)
    SO2_g_0  (Optional[float]):   initial SO2  gas volume concentration (mole fraction) [1]
    O3_g_0   (Optional[float]):   initial O3   gas volume concentration (mole fraction) [1]
    H2O2_g_0 (Optional[float]):   initial H2O2 gas volume concentration (mole fraction) [1]
    chem_sys (Optional[string]):  accepted values: 'open', 'closed'
                                  (in open/closed system gas volume concentration in the air doesn't/does change 
                                   due to chemical reactions)
    chem_dsl (Optional[bool]):    on/off for dissolving chem species into droplets
    chem_dsc (Optional[bool]):    on/off for dissociation of chem species in droplets
    chem_rct (Optional[bool]):    on/off for oxidation of S_IV to S_VI
    chem_spn (Optional[int]):     number of spinup timesteps before enabling chemical reactions
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

  micro = _micro_init(opts, state, info)
  with _output_init(micro, opts) as fout:
    # adding chem state vars
    if micro.opts_init.chem_switch:
      state.update({ "SO2_g" : SO2_g_0, "O3_g" : O3_g_0, "H2O2_g" : H2O2_g_0 })
      state.update({ "SO2_a" : 0.,      "O3_a" : 0.,     "H2O2_a" : 0.       })

    # t=0 : init & save
    _output(fout, opts, micro, state, 0)

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
      _micro_step(micro, state, info, opts, it)

      # TODO: only if user wants to stop @ RH_max
      #if (state["RH"] < info["RH_max"]): break
 
      # output
      if (it % outfreq == 0): 
        rec = it/outfreq
        _output(fout, opts, micro, state, rec)
 
    _save_attrs(fout, info)
    _save_attrs(fout, opts)

    
def _arguments_checking(args):
  if (args["gstdev"] == 1): raise Exception("standar deviation should be != 1 to avoid monodisperse distribution")
  if (args["T_0"] < 273.15): raise Exception("temperature should be larger than 0C - microphysics works only for warm clouds")
  if (args["r_0"] < 0): raise Exception("water vapour should be larger than 0")
  if (args["w"] < 0): raise Exception("vertical velocity should be larger than 0")
  if (args["kappa"] <= 0): raise Exception("kappa hygroscopicity parameter should be larger than 0 ")

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
