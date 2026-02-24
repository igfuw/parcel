#!/usr/bin/env python

import sys

from argparse import ArgumentParser, RawTextHelpFormatter
from packaging.version import Version
from scipy import __version__ as scipy_version
assert Version(scipy_version) >= Version("0.12"), "see https://github.com/scipy/scipy/pull/491"

from scipy.io import netcdf
import json, inspect, numpy as np
import pdb
import subprocess

from libcloudphxx import common, lgrngn, blk_1m
from libcloudphxx import git_revision as libcloud_version

parcel_version = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip()

# import refactored modules
from parcel_common import _Chem_g_id, _Chem_a_id, lognormal, sum_of_lognormals, _stats, _p_hydro_const_rho, _p_hydro_const_th_rv, _arguments_checking, _init_sanity_check, _w_eval
from micro_lgrngn import _micro_init as _micro_init_lgrngn, _micro_step as _micro_step_lgrngn
from micro_blk_1m import _opts_init_blk_1m, _micro_step_blk_1m
from micro_blk_1m_ice import _opts_init_blk_1m_ice, _micro_step_blk_1m_ice
from io_output import _output_bins, _output_init, _output_init_blk_1m, _output_init_blk_1m_ice, _output_save, _save_attrs, _output

def parcel(dt = .1, z_max = 200., w = 1., T_0 = 300., p_0 = 101300.,
  r_0 = -1., RH_0 = -1.,
  outfile = "test.nc",
  pprof = "pprof_piecewise_const_rhod", 
  outfreq = 100,
  scheme = "lgrngn", # microphysics scheme: lgrngn, blk_1m
  ice_switch = False,
  ice_nucl = False,
  time_dep_ice_nucl = False,
  sd_conc = 64,
  aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [60.0e6]}}',
  dry_sizes = None,
  out_bin = '{"radii": {"rght": 0.01, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "log", "left": 1e-15}}',
  SO2_g = 0., O3_g = 0., H2O2_g = 0., CO2_g = 0., HNO3_g = 0., NH3_g = 0.,
  chem_dsl = False, chem_dsc = False, chem_rct = False,
  chem_rho = 1.8e3,
  sstp_cond = 1,
  sstp_chem = 1,
  wait = 0,
  large_tail = False,
  rng_seed = None,
  rd_insol  = 0.,
  t = None,
  adaptive_sstp_cond = None,
  sstp_cond_adapt_drw2_eps = None,
  sstp_cond_adapt_drw2_max = None,
  sstp_cond_act = None,
  sstp_cond_mix = None,
  exact_sstp_cond = None,
  aerosol_independent_of_rhod = None
  ,backend = "serial"
):
  """
  Args:
    dt      (Optional[float]):    timestep [s]
    z_max   (Optional[float]):    maximum vertical displacement [m]
    t       (Optional[float|None]): simulation duration [s].
                                  Exactly one of `z_max` or `t` must be specified.
    w       (Optional[float|callable|str]): updraft velocity [m/s]
                                  - constant: number
                                  - time-dependent: callable w(t) or expression string in `t` (seconds)
                                    e.g. "1 + 0.5*np.sin(2*np.pi*t/60)"
    T_0     (Optional[float]):    initial temperature [K]
    p_0     (Optional[float]):    initial pressure [Pa]
    r_0     (Optional[float]):    initial water vapour mass mixing ratio [kg/kg]
    RH_0    (Optional[float]):    initial relative humidity
    scheme  (Optional[string]):   microphysics scheme to use: 'lgrngn', 'blk_1m'
    backend (Optional[str]):      lgrngn backend to use: 'serial', 'openmp', 'cuda' (only used when scheme='lgrngn')
    ice_switch (Optional[bool]):  enable ice microphysics
    ice_nucl (Optional[bool]):    enable ice nucleation in lagrangian scheme
    time_dep_ice_nucl (Optional[bool]): enable time-dependent ice nucleation in lagrangian scheme
    outfile (Optional[string]):   output netCDF file name
    outfreq (Optional[int]):      output interval (in number of time steps)
    pprof   (Optional[string]):   method to calculate pressure profile used to calculate
                                  dry air density that is used by the super-droplet scheme
                                  valid options are: pprof_const_th_rv, pprof_const_rhod, pprof_piecewise_const_rhod
    wait (Optional[float]):       number of timesteps to run parcel model with vertical velocity=0 at the end of simulation
                                  (added for testing)
    sd_conc (Optional[int]):      number of moving bins (super-droplets)

    
    aerosol (Optional[json str]): dict of dicts defining aerosol distribution, e.g.:

                                  {"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.02e-6, 0.07e-7], "gstdev": [1.4, 1.2], "n_tot": [120.0e6, 80.0e6]}
                                   "gccn"            : {"kappa": 1.28, "mean_r": [2e-6],             "gstdev": [1.6],      "n_tot": [1e2]}}

                                  where kappa  - hygroscopicity parameter (see doi:10.5194/acp-7-1961-2007)
                                        mean_r - lognormal distribution mean radius [m]                    (list if multimodal distribution)
                                        gstdev - lognormal distribution geometric standard deviation       (list if multimodal distribution)
                                        n_tot  - lognormal distribution total concentration under standard
                                                 conditions (T=20C, p=1013.25 hPa, rv=0) [m^-3]            (list if multimodal distribution)
                                                 
    dry_sizes (Optional[json str|dict|None]): discrete aerosol bins used to set libcloudphxx `opts_init.dry_sizes`.
                                      Can be used together with `aerosol`/dry_distros.
                                      Format example:
                                      {
                                        "ammonium_sulfate": {
                                          "kappa": 0.61,
                                          "bins": {
                                            "1e-6":  [30.0, 15],
                                            "15e-6": [10.0,  5]
                                          }
                                        }
                                      }
                                      where bins map dry_radius_m -> [STP_concentration_1_per_m3, number_of_SDs]

    large_tail (Optional[bool]) : use more SD to better represent the large tail of the initial aerosol distribution

    rng_seed (Optional[int]) :  seed for random number generator

    out_bin (Optional[json str]): dict of dicts defining spectrum diagnostics, e.g.:

                                  {"radii": {"rght": 0.0001,  "moms": [0],          "drwt": "wet", "nbin": 26, "lnli": "log", "left": 1e-09},
                                   "cloud": {"rght": 2.5e-05, "moms": [0, 1, 2, 3], "drwt": "wet", "nbin": 49, "lnli": "lin", "left": 5e-07}}
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
    chem_rho (Optional[float]):   aerosol/droplet material density for chemistry [kg/m3]
    aerosol_independent_of_rhod (Optional[bool]): on/off for initial aerosol concentration independent of rhod (assumed at STP otherwise)

    # Coalescence Substepping controls
    sstp_chem (Optional[int]):    substeps per timestep for chemistry (>=1)

    # condensation substepping controls
    sstp_cond (Optional[int]):    substeps per dynamical timestep for condensation/evaporation (>=1)
    adaptive_sstp_cond (Optional[bool]):   on/off for adaptive substepping for condensation/evaporation
    sstp_cond_adapt_drw2_eps (Optional[float]):   tolerance parameter for adaptive condensation/evaporation substepping
    sstp_cond_adapt_drw2_max (Optional[float]):   maximum relative change of rw2 for adaptive condensation/evaporation substepping
    sstp_cond_act (Optional[int]):    substeps for (de)activating droplets
    sstp_cond_mix (Optional[bool]):   on/off mixing of thermodynamic variables between superdroplets after each condensation substep
    exact_sstp_cond (Optional[bool]): on/off for per-particle condensation substepping (per-cell if off)

    # Misc
    rd_insol (Optional[float]):   insoluble dry radius offset/addition used by selected microphysics (if applicable) [m]

   """
  # packing function arguments into "opts" dictionary
  args, _, _, _ = inspect.getargvalues(inspect.currentframe())
  opts = dict()
  for k in args:
    opts[k] = locals()[k]

  # Enforce stop condition selection: either z_max or t (but not both).
  # Backwards compatible default is z_max (t=None).
  # if opts["t"] is not None and opts["z_max"] is not None:
  #   raise ValueError("Specify only one stop condition: either z_max or t (not both)")
  # if opts["t"] is None and opts["z_max"] is None:
  #   raise ValueError("You must specify exactly one stop condition: z_max or t")
  # if opts["t"] is not None and opts["t"] <= 0:
  #   raise ValueError("t must be > 0")
  # if opts["z_max"] is not None and opts["z_max"] <= 0:
  #   raise ValueError("z_max must be > 0")

  # parsing json specification of output spectra
  spectra = json.loads(opts["out_bin"])

  # parsing json specification of init aerosol spectra (if provided)
  aerosol = json.loads(opts["aerosol"]) if isinstance(opts.get("aerosol"), str) else opts.get("aerosol")

  # allow passing dry_sizes as dict or json string
  if opts.get("dry_sizes") is not None and isinstance(opts.get("dry_sizes"), str):
    dry_sizes = json.loads(opts["dry_sizes"])
  else:
    dry_sizes = opts.get("dry_sizes")

  # default water content
  if ((opts["r_0"] < 0) and (opts["RH_0"] < 0)):
    print("both r_0 and RH_0 negative, using default r_0 = 0.022")
    r_0 = .022
  # water coontent specified with RH
  if ((opts["r_0"] < 0) and (opts["RH_0"] >= 0)):
    r_0 = common.eps * opts["RH_0"] * common.p_vs(T_0) / (p_0 - opts["RH_0"] * common.p_vs(T_0))

  # sanity checks for arguments
  _arguments_checking(opts, spectra, aerosol, ice_switch)

  th_0 = T_0 * (common.p_1000 / p_0)**(common.R_d / common.c_pd)

  # Stopping condition differs for constant vs. time-dependent w.
  # For constant w: keep original behaviour (nt computed from stop condition).
  # For variable w: integrate until the stop condition is met.
  # w0 = _w_eval(w, 0.0)
  # if w0 <= 0 and isinstance(w, (int, float, np.floating)) and opts["z_max"] is not None:
  #   raise ValueError("For constant w with z_max stop, expected w>0 to reach z_max")

  # if isinstance(w, (int, float, np.floating)):
  #   if opts["t"] is not None:
  #     nt = int(np.ceil(float(opts["t"]) / dt))
  #   else:
  #     nt = int(opts["z_max"] / (float(w) * dt))
  # else:
  #   nt = None

  state = {
    "t" : 0, "z" : 0,
    "r_v" : np.array([r_0]), "p" : p_0,
    "th_d" : np.array([common.th_std2dry(th_0, r_0)]),
    "rhod" : np.array([common.rhod(p_0, th_0, r_0)]),
    "T" : None, "RH" : None
  }

  if scheme == "lgrngn" and ice_switch:
    state["ice_mix_ratio"] = np.array([0.0])

  if scheme == "blk_1m":
    state["rc"] = np.array([0.0])  # initial cloud water
    state["rr"] = np.array([0.0])  # initial rain water
    if ice_switch:
      state["ria"] = np.array([0.0])  # initial ice A
      state["rib"] = np.array([0.0])  # initial ice B

  if opts["chem_dsl"] or opts["chem_dsc"] or opts["chem_rct"]:
    for key in _Chem_g_id.keys():
      state.update({ key : np.array([opts[key]])})

  info = { "RH_max" : 0, "libcloud_Git_revision" : libcloud_version,
           "parcel_Git_revision" : parcel_version }

  _init_sanity_check(state, info)

  # Initialize the selected scheme
  if scheme == "lgrngn":
      micro = _micro_init_lgrngn(aerosol, opts, state)
      fout = _output_init(micro, opts, spectra)
  elif scheme == "blk_1m" and not ice_switch:
      micro_opts = _opts_init_blk_1m()
      fout = _output_init_blk_1m(opts)
  elif scheme == "blk_1m" and ice_switch:
      micro_opts = _opts_init_blk_1m_ice()
      fout = _output_init_blk_1m_ice(opts)
  else:
      raise ValueError("Unknown scheme type. Use 'lgrngn', 'blk_1m'.")
    
  with fout:
      # adding chem state vars - only for lgrngn scheme
      if scheme == "lgrngn" and micro.opts_init.chem_switch:
        state.update({ "SO2_a" : 0.,"O3_a" : 0.,"H2O2_a" : 0.,})
        state.update({ "CO2_a" : 0.,"HNO3_a" : 0.})

        micro.diag_all() # selecting all particles
        micro.diag_chem(_Chem_a_id["NH3_a"])
        state.update({"NH3_a": np.frombuffer(micro.outbuf())[0]})

      # t=0 : init & save
      if scheme == "lgrngn":
          _output(fout, opts, micro, state, 0, spectra)
      elif scheme == "blk_1m":
          _output_save(fout, state, 0)  # simpler output for blk_1m

      # timestepping
      rec = 0
      it = 0
      max_steps_guard = 5_000_000  # safety for pathological w(t)
      while True:
        # if nt is not None:
        #   if it >= nt:
        #     break
        # else:
        # variable-w stopping conditions
        if opts["t"] is not None:
          if state["t"] >= opts["t"]:
            break
        else:
          if state["z"] >= opts["z_max"]:
            break

        if it >= max_steps_guard:
          raise RuntimeError("Exceeded safety step limit while integrating variable w(t)")

        it += 1

        # vertical velocity at current (start-of-step) time
        w_it = _w_eval(w, state["t"])

        # diagnostics
        # the reasons to use analytic solution:
        # - independent of dt
        # - same as in 2D kinematic model
        state["z"] += w_it * dt
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
          state["p"] = _p_hydro_const_rho(w_it*dt, state["p"], state["rhod"][0])

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
        if scheme == "lgrngn":
          _micro_step_lgrngn(micro, state, info, opts)
        elif scheme == "blk_1m" and not ice_switch:
          _micro_step_blk_1m(micro_opts, state, info, opts)
        elif scheme == "blk_1m" and ice_switch:
          _micro_step_blk_1m_ice(micro_opts, state, info, opts)

        # TODO: only if user wants to stop @ RH_max
        #if (state["RH"] < info["RH_max"]): break

        # output
        if (it % outfreq == 0):
          # if nt is not None and nt > 0:
          #   print(str(round(it / (nt * 1.) * 100, 2)) + " %")
          if opts["t"] is not None:
            print(str(round(state["t"] / (opts["t"] * 1.) * 100, 2)) + " %")
          if opts["z_max"] is not None:
            print(str(round(state["z"], 1)) + " / " + str(opts["z_max"]) + " m")
          
          rec = it/outfreq
          if scheme == "lgrngn":
            _output(fout, opts, micro, state, rec, spectra)
          elif scheme == "blk_1m":
            _output_save(fout, state, rec)

      # _save_attrs(fout, info)
      _save_attrs(fout, opts)
      print("post _save_attrs")

      if wait != 0:
        for it in range (it+1, it+wait):
          state["t"] = it * dt
          if scheme == "lgrngn":
            _micro_step_lgrngn(micro, state, info, opts)
          elif scheme == "blk_1m" and not ice_switch:
            _micro_step_blk_1m(micro_opts, state, info, opts)
          elif scheme == "blk_1m" and ice_switch:
            _micro_step_blk_1m_ice(micro_opts, state, info, opts)

          if (it % outfreq == 0):
            rec = it/outfreq
            if scheme == "lgrngn":
              _output(fout, opts, micro, state, rec, spectra)
            elif scheme == "blk_1m":
              _output_save(fout, state, rec)

# ensuring that pure "import parcel" does not trigger any simulation
if __name__ == '__main__':

  # getting list of argument names and their default values
  name, _, _, dflt = inspect.getfullargspec(parcel)[0:4]
  opts = dict(list(zip(name[-len(dflt):], dflt)))

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
