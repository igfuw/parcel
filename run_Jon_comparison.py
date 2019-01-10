import os

outfile = 'plots/Jon_GCCN/Jon_GCCN.nc'
p_0     = 93844.6 
T_0     = 284.3 
#RH_0    = 0.85785275093866054 
r_0     = 0.007617065 
w       = 0.4 #0
dt      = .1
sd_conc = 100
z_max   = 600

outfreq = 150 # [s]
nt      = 0
outfreq = outfreq / dt
nt      = nt / dt

p_stp = 101325
T_stp = 273.25 + 15 


out_bin = '{"dry":    {"rght": 1e-6,  "moms": [0],       "drwt": "dry", "slct": "dry", "nbin": 500, "lnli": "log", "left": 1e-9},'+\
          '"wet":     {"rght": 40e-6, "moms": [0],       "drwt": "wet", "slct": "wet", "nbin": 500, "lnli": "log", "left": 1e-9},'+\
          '"all":     {"rght": 1e-2,  "moms": [0, 1, 2], "drwt": "wet", "slct": "wet", "nbin": 1,   "lnli": "log", "left": 1e-9},'+\
          '"cloud":   {"rght": 1e-2,  "moms": [0, 1, 2], "drwt": "wet", "slct": "wet", "nbin": 1,   "lnli": "log", "left": 1e-6},'+\
          '"rd20nm":  {"rght": 21e-9,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 19e-9},'+\
          '"rd31nm":  {"rght": 32e-9,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 30e-9},'+\
          '"rd152nm": {"rght": 153e-9,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 151e-9},'+\
          '"rd337nm": {"rght": 338e-9,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 336e-9},'+\
          '"rd500nm": {"rght": 501e-9,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 499e-9}}'

# pristine
#n_0 = 125e6
#n_1 = 65e6
#n_0_stp_corrected = n_0  * p_stp / p_0 * T_0 / T_stp
#n_1_stp_corrected = n_1  * p_stp / p_0 * T_0 / T_stp
#aerosol = '{"NaCl": {"kappa": 1.28, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": ['+str(n_0_stp_corrected)+','+str(n_1_stp_corrected)+']}}'

# modified polluted
n_0 = 48e6
n_1 = 114e6
n_0_stp_corrected = n_0  * p_stp / p_0 * T_0 / T_stp
n_1_stp_corrected = n_1  * p_stp / p_0 * T_0 / T_stp
aerosol = '{"NaCl": {"kappa": 1.28, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": ['+str(n_0_stp_corrected)+','+str(n_1_stp_corrected)+']}}'


#python parcel.py --outfile plots/Jon_GCCN/Jon_GCCN.nc --p_0 93844.6 --T_0 284.3 --RH_0 0.85785275093866054 --w 0. --dt .05 --sd_conc 1000 --z_max 788 --outfreq 10
os.system("python parcel.py --outfile "+str(outfile)+" --p_0 "+str(p_0)+" --T_0 "+str(T_0)+" --r_0 "+str(r_0)+" --w "+str(w)+" --dt "+str(dt)+" --sd_conc "+str(sd_conc)+" --z_max "+str(z_max)+" --nt "+str(int(nt))+" --aerosol \'"+str(aerosol)+"\' --out_bin \'"+str(out_bin)+"\' --outfreq "+str(int(outfreq)))
