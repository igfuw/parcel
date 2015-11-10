from parcel import parcel
from libcloudphxx import common as cm
from functions import *

"""
Common initail conditions for: test_chem_closed_dsl, test_chem_closed_dsc, test_chem_closed_rct

"""

# initial condition
RH_init = .95
T_init  = 285.2
p_init  = 95000.
r_init  = rh_to_rv(RH_init, T_init, p_init)

# calculate rhod for initial gas mixing ratio
rhod_init   = rhod_calc(T_init, p_init, r_init)
# initial condition for trace geses
SO2_g_init  = mole_frac_to_mix_ratio(200e-12, p_init, cm.M_SO2,  T_init, rhod_init)
O3_g_init   = mole_frac_to_mix_ratio(50e-9,   p_init, cm.M_O3,   T_init, rhod_init)
H2O2_g_init = mole_frac_to_mix_ratio(500e-12, p_init, cm.M_H2O2, T_init, rhod_init)
CO2_g_init  = mole_frac_to_mix_ratio(360e-6,  p_init, cm.M_CO2,  T_init, rhod_init)
NH3_g_init  = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_NH3,  T_init, rhod_init)
HNO3_g_init = mole_frac_to_mix_ratio(100e-12, p_init, cm.M_HNO3, T_init, rhod_init)

# aerosol size distribution
mean_r = .08e-6 / 2.
gstdev = 2.
n_tot  = 566.e6

# chem process toggling
chem_dsl = False
chem_dsc = False
chem_rct = False
chem_sys = 'closed'
chem_rho = 1.8e3

# output
z_max   = 600
dt      = .1
w       = .5
outfreq = int(z_max / dt / 30) 
sd_conc = 4
outfile = "TODO.nc"

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
               'T_0': T_init, 'p_0': p_init, 'r_0': r_init,\
               'SO2_g': SO2_g_init, 'O3_g':  O3_g_init, 'H2O2_g': H2O2_g_init,\
               'CO2_g': CO2_g_init, 'NH3_g': NH3_g_init,'HNO3_g': HNO3_g_init,\
               'chem_sys': chem_sys, 'chem_rho': chem_rho, 'outfile': outfile,\
               'chem_dsl': chem_dsl, 'chem_dsc': chem_dsc, 'chem_rct': chem_rct,\
               'n_tot': n_tot, 'mean_r': mean_r, 'gstdev': gstdev,\
               'sd_conc': sd_conc, 'out_bin': out_bin}
