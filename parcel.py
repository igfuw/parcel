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

Chem_ga_id = ["SO2", "H2O2", "O3"]
Chem_aq_id = Chem_ga_id + ["HSO3"]

Chem_id = {
  "SO2"  : lgrngn.chem_species_t.SO2,
  "H2O2" : lgrngn.chem_species_t.H2O2,
  "O3"   : lgrngn.chem_species_t.O3,
  "HSO3" : lgrngn.chem_species_t.HSO3                                                      
}

def micro_init(opts, state, info):
  # sanity check
  stats(state, info)
  if (state["RH"] > 1): raise Exception("Please supply initial T,p,r_v below supersaturation")

  # using nested function to get access to opts
  def lognormal(lnr):
    from math import exp, log, sqrt, pi
    return opts["n_tot"] * exp(
      -(lnr - log(opts["mean_r"]))**2 / 2 / log(opts["stdev"])**2
    ) / log(opts["stdev"]) / sqrt(2*pi);

  # lagrangian scheme options
  opts_init = lgrngn.opts_init_t()  
  for opt in ["dt", "sd_conc_mean"]:  
    setattr(opts_init, opt, opts[opt])
  opts_init.dry_distros = {opts["kappa"]:lognormal}
  opts_init.kernel = lgrngn.kernel_t.geometric #TODO: will not be needed soon (libcloud PR #89)
  opts_init.chem_switch = True 

  # initialitation
  micro = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  micro.init(state["th_d"], state["r_v"], state["rhod"])
  return micro

def micro_step(micro, state, info, chem_gas):
  libopts = lgrngn.opts_t()
  for id in chem_gas:
    libopts.chem_gas[Chem_id[id]] = chem_gas[id]

  micro.step_sync(libopts, state["th_d"], state["r_v"], state["rhod"]) 

  # new = diag_chem(Chem_id[id])
  for id in chem_gas:
    old = 0 #libopts.chem_gas[Chem_id[id]]
    
    new = 0 #np.frombuffer(micro.outbuf())    
    chem_gas[id] -= (new - old)

def stats(state, info):
  state["T"] = np.array([common.T(state["th_d"][0], state["rhod"][0])])
  state["RH"] = state["p"] * state["r_v"] / (state["r_v"] + common.eps) / common.p_vs(state["T"][0])
  info["RH_max"] = max(info["RH_max"], state["RH"])

def histo(bins, micro, opts, chem_aq):
  r_min = 0
  i = 0
  for r_max in opts["radii"]:
    micro.diag_wet_rng(r_min, r_max)

    micro.diag_wet_mom(0) # #/kg dry air
    bins["conc"][i] = np.frombuffer(micro.outbuf())

    for id in Chem_aq_id:
      micro.diag_chem(Chem_id[id])
      chem_aq[id][i] = np.frombuffer(micro.outbuf())

    r_min = r_max
    i += 1

def output_init(opts):
  # file & dimensions
  fout = netcdf.netcdf_file(opts["outfile"], 'w')
  fout.createDimension('t', None)
  fout.createDimension('radii', opts["radii"].shape[0]) #TODO: r_d, cloud only; #TODO: r_w vs. r_v - might be misleading
  
  units = {"z" : "m", "t" : "s", "r_v" : "kg/kg", "th_d" : "K", "rhod" : "kg/m3", 
    "p" : "Pa", "T" : "K", "RH" : "1", "conc" : "(kg of dry air)^-1"
  }
  for id in Chem_id:
    units[id] = "todo"

  for name, unit in units.iteritems():
    if name in Chem_aq_id + ["conc"]:
      dims = ('t','radii')
    else:
      dims = ('t',)

    fout.createVariable(name, 'd', dims)
    fout.variables[name].unit = unit

  return fout

def output_save(fout, state, rec):
  for var, val in state.iteritems():
    fout.variables[var][rec] = val

def save_attrs(fout, dictnr):
  for var, val in dictnr.iteritems():
    setattr(fout, var, val)

def output(fout, opts, micro, bins, state, chem_gas, chem_aq, rec):
  histo(bins, micro, opts, chem_aq)
  output_save(fout, state, rec)
  output_save(fout, bins, rec)
  output_save(fout, chem_aq, rec) 
  output_save(fout, chem_gas, rec)

 
def parcel(dt=.1, z_max=200, w=1, T_0=300, p_0=101300, r_0=.022, outfile="test.nc", 
  outfreq=100, sd_conc_mean=64, kappa=.5,
  mean_r = .04e-6 / 2, stdev  = 1.4, n_tot  = 60e6, 
  radii = 1e-6 * pow(10, -3 + np.arange(26) * .2), 
  SO2_0 = 44, O3_0 = 44, H2O2_0 = 44
):
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
  chem_aq = dict(zip(Chem_aq_id, len(Chem_aq_id)*[np.empty(radii.shape[0])]))
  with output_init(opts) as fout:
    # t=0 : init & save
    micro = micro_init(opts, state, info)
    output(fout, opts, micro, bins, state, chem_gas, chem_aq, 0)

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
      micro_step(micro, state, info, chem_gas)
      stats(state, info)
    
      # TODO: only if user wants to stop @ RH_max
      #if (state["RH"] < info["RH_max"]): break

      # output
      if (it % outfreq == 0): 
        rec = it/outfreq
        output(fout, opts, micro, bins, state, chem_gas, chem_aq, rec)
 
    save_attrs(fout, info)
    save_attrs(fout, opts)

# ensuring that pure "import parcel" does not trigger any simulation
if __name__ == '__main__':
  parcel()
