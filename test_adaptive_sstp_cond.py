"""
Run adaptive condensation substepping and plot results as in the MSc thesis of Piotr Bartman (Sec. 3.4, Fig. 5 therein)
"""

import sys, os
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from pathlib import Path
from typing import List

sstp_cond_max = 10
z_max = 4000.0

def run_scheme(w_max, adaptive, outfile, *, sstp_cond=sstp_cond_max):
    args = dict(
        p_0=100000,
        RH_0=0.9,
        T_0=260,
        aerosol = None,
        # aerosol = '{"pristine": {"kappa": 0.61, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}',      # aerosol=None,
        sd_conc=100,#pow(2,10),#1024,#256,
        # dry_sizes={"Bartman": {"kappa": 0.2, "bins": {
        #                                     str(r_dry):  [N_STP, 1]
        #                                   }}},
        dt=1,
        z_max=None,
        # t=300,
        w=lambda t: w_max * np.pi / 2. * np.sin(np.pi*t*w_max/z_max), # z_half = z_max
        # w=lambda t: w_max * np.pi / 2. * np.sin(np.pi*t/(z_max * w_max)), # z_half = z_max
        # r_0=0.022,
        outfile=outfile,
        # outfreq=1,
        scheme="lgrngn",
        # out_bin='{"radius": {"rght": 1, "moms": [0,1], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 1e-15}}',
        out_bin='{"cloud": {"rght": 1, "moms": [0,1], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-6}}',
        sstp_cond=sstp_cond,
        adaptive_sstp_cond=adaptive,
        # adaptive substepping parameters are injected below (only when adaptive=True)
        sstp_cond_adapt_drw2_eps=None,
        sstp_cond_adapt_drw2_max=None,
        sstp_cond_act=None,
        sstp_cond_mix   = False, #cant be True for adaptive
        exact_sstp_cond = True, # if adaptive else False,       
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = True,
        rd_insol = 0.1e-6
    )

    # NOTE: we allow passing these in through function attributes set outside.
    if hasattr(run_scheme, "aerosol"):
        args["aerosol"] = run_scheme.aerosol
    # They only make sense for adaptive_sstp_cond=True.
    if adaptive:
        if hasattr(run_scheme, "sstp_cond_adapt_drw2_eps"):
            args["sstp_cond_adapt_drw2_eps"] = float(run_scheme.sstp_cond_adapt_drw2_eps)
        if hasattr(run_scheme, "sstp_cond_adapt_drw2_max"):
            args["sstp_cond_adapt_drw2_max"] = float(run_scheme.sstp_cond_adapt_drw2_max)
        if hasattr(run_scheme, "sstp_cond_act"):
            args["sstp_cond_act"] = int(run_scheme.sstp_cond_act)        

    args["t"] = 2. * z_max / w_max # twice the time to to reach z=z_max
    args["outfreq"] = 1# args["t"] // 100  # save 100 points
    # args["t"] = 300
    # args["t"] = 1
    print("t: ", args["t"])
    parcel(**args)

    with netcdf.netcdf_file(outfile, 'r') as f:
        # rv = np.array(f.variables['r_v'][:])
        # th_d = np.array(f.variables['th_d'][:])
        z = np.array(f.variables['z'][:])
        RH = np.array(f.variables['RH'][:])
        sstp_cond_mean = np.array(f.variables['sstp_cond_mean'][:]) if 'sstp_cond_mean' in f.variables else None
        sstp_cond_mean[0] = sstp_cond_mean[1] if sstp_cond_mean is not None else None # at t=0 sstp_cond_mean=0, because its set only during the firs step (?)
        act_mom0 = np.array(f.variables['act_m0'][:]).squeeze()
        step_cond_walltime_ms = np.array(f.variables['step_cond_walltime_ms'][:]).squeeze() if 'step_cond_walltime_ms' in f.variables else None
    return RH, z, sstp_cond_mean, act_mom0, step_cond_walltime_ms

# --- batch scenarios ---

# baseline - basically no adaptation, very relaxed conditions
baseline = dict(
    eps=1e6, #1e-1,
    max=1e6, #100,
    act=1,  # 1 means disabled
)

vary_eps = [1e-1, 1e-2, 1e-3]

def make_figure(aerosol_name, aerosol, xmax):
    run_scheme.aerosol = aerosol
    # rows: w_max; cols: eps
    # w_max_list = [0.1, 1., 2.5, 5.0]
    w_max_list = [5.0]
    fig, axes = plt.subplots(len(w_max_list), len(vary_eps), figsize=(15.0, 15.0), sharex=True, sharey=True, squeeze=False)

    generated_nc_files: List[str] = []

    # shared colormap settings for sstp_cond_dt (= sstp_cond_mean here)
    cmap_dt = "gnuplot"
    norm_dt = mcolors.Normalize(vmin=1, vmax=sstp_cond_max)

    for i, w_max in enumerate(w_max_list):
        # --- reference run (non-adaptive) once per w_max ---
        outfile_ref = f"test_adaptive_sstp_cond_{aerosol_name}_w{w_max:g}_ref_adapt0.nc"
        RH_ref, z_ref, _, act_mom0_ref, step_cond_ref_ms = run_scheme(w_max, False, outfile_ref)
        generated_nc_files.append(outfile_ref)
        x_ref = act_mom0_ref / 1e6
        y_ref = z_ref

        ref_step_cond_mean_ms = float(np.nanmean(step_cond_ref_ms)) if step_cond_ref_ms is not None else float("nan")

        for j, eps in enumerate(vary_eps):
            ax = axes[i, j]

            # overlay reference
            ax.plot(x_ref, y_ref, color="0.6", linewidth=2.0, zorder=1)

            run_scheme.sstp_cond_adapt_drw2_eps = eps
            run_scheme.sstp_cond_adapt_drw2_max = baseline["max"]
            run_scheme.sstp_cond_act = baseline["act"]

            outfile = f"test_adaptive_sstp_cond_{aerosol_name}_w{w_max:g}_eps{eps:.0e}_adapt1.nc"
            RH, z, sstp_cond_mean, act_mom0, step_cond_walltime_ms = run_scheme(w_max, True, outfile)
            generated_nc_files.append(outfile)

            x = act_mom0 / 1e6
            y = z

            sstp_cond_avg = float(np.nanmean(sstp_cond_mean)) if sstp_cond_mean is not None else float("nan")
            step_cond_mean_ms = float(np.nanmean(step_cond_walltime_ms)) if step_cond_walltime_ms is not None else float("nan")
            speedup = (ref_step_cond_mean_ms / step_cond_mean_ms) if (np.isfinite(ref_step_cond_mean_ms) and np.isfinite(step_cond_mean_ms) and step_cond_mean_ms > 0) else float("nan")

            ax.text(
                0.98,
                0.05,
                f"sstp_cond_avg={sstp_cond_avg:.2f}\nref={ref_step_cond_mean_ms:.2f} ms\nspeedup={speedup:.2f}x",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7),
                zorder=5,
            )

            if sstp_cond_mean is None:
                ax.plot(x, y, linewidth=2.0, color="C0", zorder=2)
            else:
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments, cmap=cmap_dt, norm=norm_dt)
                lc.set_array(sstp_cond_mean[1:])
                lc.set_linewidth(2.0)
                lc.set_zorder=2
                ax.add_collection(lc)

            if i == 0:
                ax.set_title(f"eps={eps:.0e}")
            if j == 0:
                ax.set_ylabel(f"w_max={w_max:g}\nHeight [m]")

    for i in range(len(w_max_list)):
        for j in range(len(vary_eps)):
            ax = axes[i, j]
            ax.set_xlim(-5, xmax)
            ax.set_ylim(0, z_max + 25)
            if i == len(w_max_list) - 1:
                ax.set_xlabel("activated droplets [1/mg]")

    fig.suptitle("Adaptive substepping, " + aerosol_name)
    fig.tight_layout(rect=(0, 0.10, 1, 0.97))

    # shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_dt, norm=norm_dt)
    sm.set_array([])
    cax = fig.add_axes([0.15, 0.04, 0.70, 0.025])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("sstp_cond_mean [1]")

    out_png = "test_adaptive_sstp_cond_"+aerosol_name+".png"
    plt.savefig(out_png, dpi=200)

    # cleanup generated NetCDF files for this figure
    # for p in generated_nc_files:
    #     try:
    #         Path(p).unlink(missing_ok=True)
    #     except TypeError:
    #         try:
    #             if Path(p).exists():
    #                 Path(p).unlink()
    #         except OSError:
    #             pass
    #     except OSError:
    #         pass

    return fig

make_figure('pristine', '{"DYCOMS": {"kappa": 0.61, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}', 200)
make_figure('polluted', '{"polluted": {"kappa": 0.61, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}', 600)
# plt.show()
