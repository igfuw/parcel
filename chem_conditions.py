from parcel import parcel
from libcloudphxx import common as cm
import functions as fn

"""
Common initail conditions for: 
test_chem_closed_dsl, test_chem_closed_dsc,    test_chem_closed_rct
tes_chem_mass_const,  test_chem_init_spectrum, test_sd_convergence

"""
# initial conditions (from Kreidenweis et al 2003)
RH_init = .95
T_init  = 285.2
p_init  = 95000.
r_init  = fn.rh_to_rv(RH_init, T_init, p_init)

# STP conditions (needed to initialize dry radii distribution)
p_stp = 101325
T_stp = 273.15 + 15

# calculate rhod for initial gas mixing ratio
rhod_init   = fn.rhod_calc(T_init, p_init, r_init)

# calculate density of air in the model and in standard conditions 
# (needed to initialize dry radii distribution)
rho_init = p_init / T_init / (r_init / (1.+r_init) * cm.R_v + 1./ (1.+r_init) * cm.R_d)
rho_stp  = p_stp  / T_stp / cm.R_d

# initial condition for trace geses
SO2_g_init  = fn.mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
O3_g_init   = fn.mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
H2O2_g_init = fn.mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
CO2_g_init  = fn.mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
NH3_g_init  = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
HNO3_g_init = fn.mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

# aerosol size distribution
mean_r = .04e-6
gstdev = 2.
n_tot  = 566.e6 * rho_stp / rho_init

# chem process toggling
chem_dsl = False
chem_dsc = False
chem_rct = False
chem_sys = 'closed'
chem_rho = 1.8e3

# output
z_max   = 1400
dt      = 1
w       = .5
outfreq = int(z_max / dt / 30) 
sd_conc = 4
outfile = "TODO.nc"

#substeps for condensation and chemistry
sstp_cond = 10

# output moments for the chemistry quicklook plot (the same for all tests)
out_bin = '{"plt_rw":   {"rght": 1,    "left":    0, "drwt": "wet", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
            "plt_rd":   {"rght": 1,    "left":    0, "drwt": "dry", "lnli": "lin", "nbin": 1, "moms": [0, 1, 3]},\
            "plt_ch":   {"rght": 1,    "left":    0, "drwt": "dry", "lnli": "lin", "nbin": 1,\
                         "moms": ["O3_a",   "H2O2_a", "H", "OH",\
                                 "SO2_a",  "HSO3_a", "SO3_a", "HSO4_a", "SO4_a",  "S_VI",\
                                 "CO2_a",  "HCO3_a", "CO3_a",\
                                 "NH3_a",  "NH4_a",  "HNO3_a", "NO3_a"]}}'

# saving parcel options as a dictionary
parcel_dict = {'dt': dt, 'z_max': z_max, 'outfreq': outfreq, 'w': w,\
               'sstp_cond': sstp_cond,\
               'T_0': T_init, 'p_0': p_init, 'r_0': r_init,\
               'SO2_g': SO2_g_init, 'O3_g':  O3_g_init, 'H2O2_g': H2O2_g_init,\
               'CO2_g': CO2_g_init, 'NH3_g': NH3_g_init,'HNO3_g': HNO3_g_init,\
               'chem_sys': chem_sys, 'chem_rho': chem_rho, 'outfile': outfile,\
               'chem_dsl': chem_dsl, 'chem_dsc': chem_dsc, 'chem_rct': chem_rct,\
               'n_tot': n_tot, 'mean_r': mean_r, 'gstdev': gstdev,\
               'sd_conc': sd_conc, 'out_bin': out_bin}
