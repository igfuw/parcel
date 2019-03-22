#!/usr/bin/env python

# TEMP TODO TEMP TODO !!!
import sys
#sys.path.insert(0, "../libcloudphxx/build/bindings/python/")
#sys.path.insert(0, "../../../libcloudphxx/build/bindings/python/")
#sys.path.insert(0, "/usr/local/lib/site-python/")
#sys.path.insert(0, "/mnt/local/pdziekan/usr/local/lib")
#sys.path.insert(0, "/mnt/local/pdziekan/lib/python/site-packages/")
sys.path.insert(0, "/usr/local/lib/python2.7/dist-packages/")
# TEMP TODO TEMP TODO !!!

from argparse import ArgumentParser, RawTextHelpFormatter

from distutils.version import StrictVersion
from scipy import __version__ as scipy_version
assert StrictVersion(scipy_version) >= StrictVersion("0.12"), "see https://github.com/scipy/scipy/pull/491"

from scipy.io import netcdf
import json, inspect, numpy as np
import pdb
import subprocess

from libcloudphxx import common, lgrngn
from libcloudphxx import git_revision as libcloud_version

parcel_version = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip()

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

def _micro_init(aerosol, aerosol_sizes, opts, state, info):

  # lagrangian scheme options
  opts_init = lgrngn.opts_init_t()
  for opt in ["dt", "sd_conc", "chem_rho", "sstp_cond"]:
    setattr(opts_init, opt, opts[opt])
  #opts_init.n_sd_max = int(1.1*opts_init.sd_conc)
  opts_init.n_sd_max = opts_init.sd_conc + 38  # some more space for the dry_sizes GCCNs

  opts_init.kernel = lgrngn.kernel_t.hall_davis_no_waals
  opts_init.terminal_velocity = lgrngn.vt_t.khvorostyanov_nonspherical;

  # read in the initial aerosol size distribution
  dry_distros = {}
  for name, dct in aerosol.iteritems(): # loop over kappas
    lognormals = []
    for i in range(len(dct["mean_r"])):
      lognormals.append(lognormal(dct["mean_r"][i], dct["gstdev"][i], dct["n_tot"][i]))
    dry_distros[dct["kappa"]] = sum_of_lognormals(lognormals)
  opts_init.dry_distros = dry_distros

  # read in the initial aerosol dry sizes
  if (aerosol_sizes != None):
    print aerosol_sizes
    dry_sizes = {}
    for name, dct in aerosol_sizes.iteritems(): # loop over kappas
      print dct["kappa"]
      print dct["r"]
      print dct["n_tot"]
      print dct["multi"]
      sizes = {}
      for i in range(len(dct["r"])):
        sizes[dct["r"][i]] = [dct["n_tot"][i], dct["multi"][i]]
      dry_sizes[dct["kappa"]] = sizes
    opts_init.dry_sizes = dry_sizes

  # better resolution for the SD tail
  if opts["large_tail"]:
      opts_init.sd_conc_large_tail = 1
      opts_init.n_sd_max = int(1e6)  # some more space for the tail SDs

  # switch off sedimentation and collisions
  opts_init.sedi_switch = False
  opts_init.coal_switch = True

  # switching on chemistry if either dissolving, dissociation or reactions are chosen
  opts_init.chem_switch = False
  if opts["chem_dsl"] or opts["chem_dsc"] or opts["chem_rct"]:
    opts_init.chem_switch = True
    opts_init.sstp_chem = opts["sstp_chem"]

  # initialisation
  micro = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  ambient_chem = {}
  if micro.opts_init.chem_switch:
    ambient_chem = dict((v, state[k]) for k,v in _Chem_g_id.iteritems())
  micro.init(state["th_d"], state["r_v"], state["rhod"], ambient_chem=ambient_chem)

  # sanity check
  _stats(state, info, micro)
  if (state["RH"] > 1): raise Exception("Please supply initial T,p,r_v below supersaturation")

  return micro

def _micro_step(micro, state, info, opts, it, fout):
  libopts = lgrngn.opts_t()
  libopts.cond = True
  libopts.coal = True
  libopts.adve = False
  libopts.sedi = False

  # chemical options
  if micro.opts_init.chem_switch:
    # chem processes: dissolving, dissociation, reactions
    libopts.chem_dsl = opts["chem_dsl"]
    libopts.chem_dsc = opts["chem_dsc"]
    libopts.chem_rct = opts["chem_rct"]

  # get trace gases
  ambient_chem = {}
  if micro.opts_init.chem_switch:
    ambient_chem = dict((v, state[k]) for k,v in _Chem_g_id.iteritems())

  # call libcloudphxx microphysics
  micro.step_sync(libopts, state["th_d"], state["r_v"], state["rhod"], ambient_chem=ambient_chem)
  micro.step_async(libopts)

  # update state after microphysics (needed for below update for chemistry)
  _stats(state, info, micro)

  # update in state for aqueous chem (TODO do we still want to have aq chem in state?)
  if micro.opts_init.chem_switch:
    micro.diag_all() # selecting all particles
    for id_str, id_int in _Chem_g_id.iteritems():
      # save changes due to chemistry
      micro.diag_chem(id_int)
      state[id_str.replace('_g', '_a')] = np.frombuffer(micro.outbuf())[0]

def _stats(state, info, micro):
  state["T"] = np.array([common.T(state["th_d"][0], state["rhod"][0])])
  state["RH"] = state["p"] * state["r_v"] / (state["r_v"] + common.eps) / common.p_vs(state["T"][0])
  micro.diag_RH()
  state["RH_libcloud"] =  np.frombuffer(micro.outbuf())[0]
  info["RH_max"] = max(info["RH_max"], state["RH"])

def _output_bins(fout, t, micro, opts, spectra):
  for dim, dct in spectra.iteritems():
    for bin in range(dct["nbin"]):
      if dct["slct"] == 'wet':
	micro.diag_wet_rng(
	  fout.variables[dim+"_r_wet"][bin],
	  fout.variables[dim+"_r_wet"][bin] + fout.variables[dim+"_dr_wet"][bin]
	)
      elif dct["slct"] == 'dry':
	micro.diag_dry_rng(
	  fout.variables[dim+"_r_dry"][bin],
	  fout.variables[dim+"_r_dry"][bin] + fout.variables[dim+"_dr_dry"][bin]
	)
      else: raise Exception("slct should be wet or dry")

      for vm in dct["moms"]:
        if type(vm) == int:
          # calculating moments
          if dct["drwt"] == 'wet':
            micro.diag_wet_mom(vm)
          elif dct["drwt"] == 'dry':
            micro.diag_dry_mom(vm)
          else: raise Exception("drwt should be wet or dry")
          fout.variables[dim+'_'+dct["drwt"]+'_mom'+str(vm)][int(t), int(bin)] = np.frombuffer(micro.outbuf())
        else:
          # calculate chemistry
          micro.diag_chem(_Chem_a_id[vm])
          fout.variables[dim+'_'+vm][int(t), int(bin)] = np.frombuffer(micro.outbuf())

def _output_init(micro, opts, spectra):
  # file & dimensions
  fout = netcdf.netcdf_file(opts["outfile"], 'w')
  fout.createDimension('t', None)
  for name, dct in spectra.iteritems():
    fout.createDimension(name, dct["nbin"])

    tmp = name + '_r_' + dct["slct"]
    fout.createVariable(tmp, 'd', (name,))
    fout.variables[tmp].unit = "m"
    fout.variables[tmp].description = "particle "+dct["slct"]+" radius (left bin edge)"

    tmp = name + '_dr_' + dct["slct"]
    fout.createVariable(tmp, 'd', (name,))
    fout.variables[tmp].unit = "m"
    fout.variables[tmp].description = "bin width"

    if dct["lnli"] == 'log':
      from math import exp, log
      dlnr = (log(dct["rght"]) - log(dct["left"])) / dct["nbin"]
      allbins = np.exp(log(dct["left"]) + np.arange(dct["nbin"]+1) * dlnr)
      fout.variables[name+'_r_'+dct["slct"]][:] = allbins[0:-1]
      fout.variables[name+'_dr_'+dct["slct"]][:] = allbins[1:] - allbins[0:-1]
    elif dct["lnli"] == 'lin':
      dr = (dct["rght"] - dct["left"]) / dct["nbin"]
      fout.variables[name+'_r_'+dct["slct"]][:] = dct["left"] + np.arange(dct["nbin"]) * dr
      fout.variables[name+'_dr_'+dct["slct"]][:] = dr
    else: raise Exception("lnli should be log or lin")

    for vm in dct["moms"]:
      if (vm in _Chem_a_id):
      	fout.createVariable(name+'_'+vm, 'd', ('t',name))
      	fout.variables[name+'_'+vm].unit = 'kg of chem species dissolved in cloud droplets (kg of dry air)^-1'
      else:
        assert(type(vm)==int)
	fout.createVariable(name+'_'+dct["drwt"]+'_mom'+str(vm), 'd', ('t',name))
	fout.variables[name+'_'+dct["drwt"]+'_mom'+str(vm)].unit = 'm^'+str(vm)+' (kg of dry air)^-1'

  units = {"z"  : "m",     "t"   : "s",     "r_v"  : "kg/kg", "th_d" : "K", "rhod" : "kg/m3",
           "p"  : "Pa",    "T"   : "K",     "RH"   : "1"    , "RH_libcloud" : "1", "curr_w" : "m/s"
  }

  if micro.opts_init.chem_switch:
    for id_str in _Chem_g_id.iterkeys():
      units[id_str] = "gas mixing ratio [kg / kg dry air]"
      units[id_str.replace('_g', '_a')] = "kg of chem species (both undissociated and ions) dissolved in cloud droplets (kg of dry air)^-1"

  for var_name, unit in units.iteritems():
    fout.createVariable(var_name, 'd', ('t',))
    fout.variables[var_name].unit = unit

  return fout

def _output_save(fout, state, rec):
  for var, val in state.iteritems():
    fout.variables[var][int(rec)] = val

def _save_attrs(fout, dictnr):
  for var, val in dictnr.iteritems():
    setattr(fout, var, val)

def _output(fout, opts, micro, state, rec, spectra):
  _output_bins(fout, rec, micro, opts, spectra)
  _output_save(fout, state, rec)

def _p_hydro_const_rho(dz, p, rho):
  # hydrostatic pressure assuming constatnt density
  return p - rho * common.g * dz

def _p_hydro_const_th_rv(z_lev, p_0, th_std, r_v, z_0=0.):
  # hydrostatic pressure assuming constatnt theta and r_v
  return common.p_hydro(z_lev, th_std, r_v, z_0, p_0)

def parcel(dt=.1, nt=0, z_max=200., z_min=-1., w=1., T_0=300., p_0=101300.,
  r_0=-1., RH_0=-1., #if none specified, the default will be r_0=.022,
  outfile="test.nc",
  pprof="pprof_piecewise_const_rhod",
  outfreq=100, sd_conc=64,
  aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [60.0e6]}}',
  aerosol_sizes = 'null', 
  out_bin = '{"radii": {"rght": 0.0001, "moms": [0], "drwt": "wet", "slct": "wet", "nbin": 1, "lnli": "log", "left": 1e-09}}',
  SO2_g = 0., O3_g = 0., H2O2_g = 0., CO2_g = 0., HNO3_g = 0., NH3_g = 0.,
  chem_dsl = False, chem_dsc = False, chem_rct = False,
  chem_rho = 1.8e3,
  sstp_cond = 1,
  sstp_chem = 1,
  wait = 0,
  top_stop = 0.,
  large_tail = False
):
  """
  Args:
    dt      (Optional[float]):    timestep [s]
    nt      (Optional[int]):      number of timesteps, only for updraft
    z_max   (Optional[float]):    maximum vertical displacement [m]
    z_min   (Optional[float]):    if non-negative, parcel will descend to z_min after reaching z_max [m]
    w       (Optional[float]):    updraft velocity [m/s]
    T_0     (Optional[float]):    initial temperature [K]
    p_0     (Optional[float]):    initial pressure [Pa]
    r_0     (Optional[float]):    initial water vapour mass mixing ratio [kg/kg]
    RH_0    (Optional[float]):    initial relative humidity
    outfile (Optional[string]):   output netCDF file name
    outfreq (Optional[int]):      output interval (in number of time steps)
    pprof   (Optional[string]):   method to calculate pressure profile used to calculate
                                  dry air density that is used by the super-droplet scheme
                                  valid options are: pprof_const_th_rv, pprof_const_rhod, pprof_piecewise_const_rhod
    wait (Optional[float]):       number of timesteps to run parcel model with vertical velocity=0 at the end of simulation
                                  (added for testing)
    top_stop (Optional[float]):   time the parcel waits at the top of the cloud before starting the descent [s]
    sd_conc (Optional[int]):      number of moving bins (super-droplets)

    aerosol (Optional[json str]): dict of dicts defining aerosol distribution, e.g.:

                                  {"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.02e-6, 0.07e-7], "gstdev": [1.4, 1.2], "n_tot": [120.0e6, 80.0e6]}
                                   "gccn"            : {"kappa": 1.28, "mean_r": [2e-6],             "gstdev": [1.6],      "n_tot": [1e2]}}

                                  where kappa  - hygroscopicity parameter (see doi:10.5194/acp-7-1961-2007)
                                        mean_r - lognormal distribution mean radius [m]                    (list if multimodal distribution)
                                        gstdev - lognormal distribution geometric standard deviation       (list if multimodal distribution)
                                        n_tot  - lognormal distribution total concentration under standard
                                                 conditions (T=20C, p=1013.25 hPa, rv=0) [m^-3]            (list if multimodal distribution)

    aerosol_sizes (Optional[json str]): dict of dicts defining aerosol by size-concentration pair, e.g.:

                                  {gccn"            : {"kappa": 1.28, "r": [2e-6, 4e-6],      "n_tot": [1e2, 5e1]}}

                                  where kappa  - hygroscopicity parameter (see doi:10.5194/acp-7-1961-2007)
                                        r -      dry radius [m]                    (list if multimodal distribution)
                                        n_tot  - total concentration under standard
                                                 conditions (T=20C, p=1013.25 hPa, rv=0) [m^-3]            (list if multimodal distribution)

    large_tail (Optional[bool]) : use more SD to better represent the large tail of the initial aerosol distribution

    out_bin (Optional[json str]): dict of dicts defining spectrum diagnostics, e.g.:

                                  {"radii": {"rght": 0.0001,  "moms": [0],          "drwt": "wet", "slct": "wet", "nbin": 26, "lnli": "log", "left": 1e-09},
                                   "cloud": {"rght": 2.5e-05, "moms": [0, 1, 2, 3], "drwt": "wet", "slct": "wet", "nbin": 49, "lnli": "lin", "left": 5e-07}}
                                  will generate five output spectra:
                                  - 0-th spectrum moment for 26 bins spaced logarithmically between 0 and 1e-4 m for dry radius
                                  - 0,1,2 & 3-rd moments for 49 bins spaced linearly between .5e-6 and 25e-6 for wet radius

                                  It can also define spectrum diagnostics for chemical compounds, e.g.:

                                  {"chem" : {"rght": 1e-6, "left": 1e-10, "drwt": "dry", "lnli": "log", "nbin": 100, "moms": ["S_VI", "NH4_a"]}}
                                  will output the total mass of H2SO4  and NH4 ions in each sizedistribution bin

                                  Valid "moms" for chemistry are:
                                    "O3_a",  "H2O2_a", "H",
                                    "SO2_a",  "S_VI",
                                    "CO2_a",
                                    "NH3_a", "HNO3_a",

    SO2_g    (Optional[float]):   initial SO2  gas mixing ratio [kg / kg dry air]
    O3_g     (Optional[float]):   initial O3   gas mixing ratio [kg / kg dry air]
    H2O2_g   (Optional[float]):   initial H2O2 gas mixing ratio [kg / kg dry air]
    CO2_g    (Optional[float]):   initial CO2  gas mixing ratio [kg / kg dry air]
    NH3_g     (Optional[float]):  initial NH3  gas mixing ratio [kg / kg dry air]
    HNO3_g   (Optional[float]):   initial HNO3 gas mixing ratio [kg / kg dry air]
    chem_dsl (Optional[bool]):    on/off for dissolving chem species into droplets
    chem_dsc (Optional[bool]):    on/off for dissociation of chem species in droplets
    chem_rct (Optional[bool]):    on/off for oxidation of S_IV to S_VI

}


   """
  # packing function arguments into "opts" dictionary
  args, _, _, _ = inspect.getargvalues(inspect.currentframe())
  opts = dict(zip(args, [locals()[k] for k in args]))

  # parsing json specification of output spectra
  spectra = json.loads(opts["out_bin"])

  # parsing json specification of init aerosol spectra
  aerosol = json.loads(opts["aerosol"])
  aerosol_sizes = json.loads(opts["aerosol_sizes"])

  # default water content
  if ((opts["r_0"] < 0) and (opts["RH_0"] < 0)):
    print "both r_0 and RH_0 negative, using default r_0 = 0.022"
    r_0 = .022
  # water coontent specified with RH
  if ((opts["r_0"] < 0) and (opts["RH_0"] >= 0)):
    r_0 = common.eps * opts["RH_0"] * common.p_vs(T_0) / (p_0 - opts["RH_0"] * common.p_vs(T_0))

  # sanity checks for arguments
  _arguments_checking(opts, spectra, aerosol)

  th_0 = T_0 * (common.p_1000 / p_0)**(common.R_d / common.c_pd)
  if(z_max > 0):
    nt_ud = int(round(z_max / (w * dt)))
    nt_dd = 0
    if(z_min >= 0):
      nt_dd = int(round((z_max - z_min) / (w * dt)))
    nt_topstop = int(round(top_stop / dt))
    nt = nt_ud + nt_topstop + nt_dd


  state = {
    "t" : 0, "z" : 0,
    "r_v" : np.array([r_0]), "p" : p_0,
    "th_d" : np.array([common.th_std2dry(th_0, r_0)]),
    "rhod" : np.array([common.rhod(p_0, th_0, r_0)]),
    "T" : None, "RH" : None, "RH_libcloud" : None,
    "curr_w" : w
  }

  if opts["chem_dsl"] or opts["chem_dsc"] or opts["chem_rct"]:
    for key in _Chem_g_id.iterkeys():
      state.update({ key : np.array([opts[key]])})

  info = { "RH_max" : 0, "libcloud_Git_revision" : libcloud_version,
           "parcel_Git_revision" : parcel_version }

  micro = _micro_init(aerosol, aerosol_sizes, opts, state, info)

  with _output_init(micro, opts, spectra) as fout:
    # adding chem state vars
    if micro.opts_init.chem_switch:
      state.update({ "SO2_a" : 0.,"O3_a" : 0.,"H2O2_a" : 0.,})
      state.update({ "CO2_a" : 0.,"HNO3_a" : 0.})

      micro.diag_all() # selecting all particles
      micro.diag_chem(_Chem_a_id["NH3_a"])
      state.update({"NH3_a": np.frombuffer(micro.outbuf())[0]})

    # t=0 : init & save
    _output(fout, opts, micro, state, 0, spectra)

    # timestepping
    for it in range(1,nt+1):
      # diagnostics
      # the reasons to use analytic solution:
      # - independent of dt
      # - same as in 2D kinematic model
     
      state["curr_w"] = w # for runs with nt specified and z_max = 0
      if(it < nt_ud):
        state["curr_w"] = w
      elif(it < nt_ud + nt_topstop):
        state["curr_w"] = 0.
      else:
        state["curr_w"] = -w

      state["z"] += state["curr_w"] * dt
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
        state["p"] = _p_hydro_const_rho(state["curr_w"]*dt, state["p"], state["rhod"][0])

      else: raise Exception("pprof should be pprof_const_th_rv, pprof_const_rhod, or pprof_piecewise_const_rhod")

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
      _micro_step(micro, state, info, opts, it, fout)

      # TODO: only if user wants to stop @ RH_max
      #if (state["RH"] < info["RH_max"]): break

      # output
      if (it % outfreq == 0):
        print str(round(it / (nt * 1.) * 100, 2)) + " %"
        rec = it/outfreq
        _output(fout, opts, micro, state, rec, spectra)

    _save_attrs(fout, info)
    _save_attrs(fout, opts)

    if wait != 0:
      for it in range (nt+1, nt+wait):
        state["t"] = it * dt
        _micro_step(micro, state, info, opts, it, fout)

        if (it % outfreq == 0):
          rec = it/outfreq
          _output(fout, opts, micro, state, rec, spectra)

def _arguments_checking(opts, spectra, aerosol):
  if opts["T_0"] < 273.15:
    raise Exception("temperature should be larger than 0C - microphysics works only for warm clouds")
  elif ((opts["r_0"] >= 0) and (opts["RH_0"] >= 0)):
    raise Exception("both r_0 and RH_0 specified, please use only one")
  elif ((opts["z_max"] > 0) and (opts["nt"] > 0)):
    raise Exception("both z_max and nt specified, please use only one")
  elif ((opts["z_min"] > 0) and (opts["nt"] > 0)):
    raise Exception("both z_min and nt specified, please use only one")
  elif ((opts["top_stop"] > 0) and (opts["nt"] > 0)):
    raise Exception("both top_stop and nt specified, please use only one")
  elif ((opts["top_stop"] < 0)):
    raise Exception("top_stop < 0")
  elif ((opts["z_max"] <= 0) and (opts["nt"] <= 0)):
    raise Exception("either z_max or nt should be larger than 0")
  elif ((opts["z_max"] > 0) and (opts["z_min"] >= opts["z_max"])):
    raise Exception("z_min >= z_max !")
  if opts["w"] < 0:
    raise Exception("vertical velocity should be larger than 0")

  for name, dct in aerosol.iteritems():
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

  for name, dct in spectra.iteritems():
    # TODO: check if name is valid netCDF identifier
    # (http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM/Identifiers.html)
    keys = ["left", "rght", "nbin", "drwt", "slct", "lnli", "moms"]
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
    if dct["slct"] not in ["dry", "wet"]:
        raise Exception(">>slct<< key in out_bin["+ name +"] must be either >>dry<< or >>wet<<")
    if dct["lnli"] not in ["lin", "log"]:
        raise Exception(">>lnli<< key in out_bin["+ name +"] must be either >>lin<< or >>log<<")
    if type(dct["nbin"]) != int:
        raise Exception(">>nbin<< key in out_bin["+ name +"] must be an integer number")
    if type(dct["moms"]) != list:
        raise Exception(">>moms<< key in out_bin["+ name +"] must be a list")
    for mom in dct["moms"]:
        if (type(mom) != int):
          if (mom not in _Chem_a_id.keys()):
            raise Exception(">>moms<< key in out_bin["+ name +"] must be a list of integer numbers or valid chemical compounds (" +str(_Chem_a_id.keys()) + ")")

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
