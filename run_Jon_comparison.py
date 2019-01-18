import os

outfile = 'plots/Jon_GCCN/Jon_GCCN.nc'
p_0     = 93844.6 
T_0     = 284.3 
#RH_0    = 0.85785275093866054 
r_0     = 0.007617065 
w       = 0.4 #0
dt      = .1
sd_conc = 10000
z_max   = 600
z_min   = 300
top_stop = 0 #[s]

outfreq = 1 # [s]
nt      = 0
outfreq = outfreq / dt
nt      = nt / dt

p_stp = 101325
T_stp = 273.25 + 15 


#out_bin = '{"dry":    {"rght": 1e-6,  "moms": [0],       "drwt": "dry", "slct": "dry", "nbin": 500, "lnli": "log", "left": 1e-9},'+\
#          '"wet":     {"rght": 40e-6, "moms": [0],       "drwt": "wet", "slct": "wet", "nbin": 500, "lnli": "log", "left": 1e-9},'+\
out_bin = '{"all":    {"rght": 1e-2,  "moms": [0, 1, 2], "drwt": "wet", "slct": "wet", "nbin": 1,   "lnli": "log", "left": 1e-9},'+\
          '"cloud":   {"rght": 1e-2,  "moms": [0, 1, 2], "drwt": "wet", "slct": "wet", "nbin": 1,   "lnli": "log", "left": 1e-6},'+\
          '"rd1.0um": {"rght": 1.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 0.99e-6},'+\
          '"rd2.0um": {"rght": 2.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 1.99e-6},'+\
          '"rd3.0um": {"rght": 3.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 2.99e-6},'+\
          '"rd4.0um": {"rght": 4.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 3.99e-6},'+\
          '"rd5.0um": {"rght": 5.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 4.99e-6},'+\
          '"rd6.0um": {"rght": 6.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 5.99e-6},'+\
          '"rd7.0um": {"rght": 7.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 6.99e-6},'+\
          '"rd8.0um": {"rght": 8.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 7.99e-6},'+\
          '"rd9.0um": {"rght": 9.01e-6,  "moms": [0, 1],    "drwt": "wet", "slct": "dry", "nbin": 1, "lnli": "log", "left": 8.99e-6},'+\
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

# GCCNs
sizes = [0.8,\
1.0, 1.2, 1.4, 1.6, 1.8,\
2.0, 2.2, 2.4, 2.6, 2.8,\
3.0, 3.2, 3.4, 3.6, 3.8,\
4.0, 4.2, 4.4, 4.6, 4.8,\
5.0, 5.2, 5.4, 5.6, 5.8,\
6.0, 6.2, 6.4, 6.6, 6.8,\
7.0, 7.2, 7.4, 7.6, 7.8,\
8.6, 9.0]

concentrations = [ 111800.,\
68490., 38400., 21820., 13300., 8496.,\
 5486., 3805. , 2593. , 1919. , 1278.,\
 988.4, 777.9 , 519.5 , 400.5 , 376.9,\
 265.3, 212.4 , 137.8 , 121.4 , 100.9,\
 122.2, 50.64 , 38.30 , 55.47 , 21.45,\
 12.95, 43.23 , 26.26 , 30.50 , 4.385,\
 4.372, 4.465 , 4.395 , 4.427 , 4.411,\
 4.522, 4.542]

multiplicities = [ 11180.,\
6849., 3840., 2182., 1330., 849.,\
 548., 380. , 259. , 191. , 127.,\
 98, 77 , 51 , 40 , 37,\
 26, 21 , 13 , 12 , 10,\
 12, 5 , 3 , 5 , 2,\
 1, 4 , 2 , 3 , 1,\
 1, 1 , 1 , 1 , 1,\
 1, 1]

aerosol_sizes = '{"NaCl": {"kappa": 1.28, "r": ['
for size in sizes:
  aerosol_sizes += str(size*1e-6)+','
#remove last ,
aerosol_sizes = aerosol_sizes[:-1]

aerosol_sizes += '], "n_tot": ['
for conc in concentrations:
  aerosol_sizes += str(conc * p_stp / p_0 * T_0 / T_stp)+','
#remove last ,
aerosol_sizes = aerosol_sizes[:-1]

aerosol_sizes += '], "multi": ['
#for multi in multiplicities:
#  aerosol_sizes += str(int(multi))+','
for conc in concentrations:
  aerosol_sizes += str(int(conc * p_stp / p_0 * T_0 / T_stp))+','
#remove last ,
aerosol_sizes = aerosol_sizes[:-1]

aerosol_sizes += ']}}'

print aerosol_sizes


#python parcel.py --outfile plots/Jon_GCCN/Jon_GCCN.nc --p_0 93844.6 --T_0 284.3 --RH_0 0.85785275093866054 --w 0. --dt .05 --sd_conc 1000 --z_max 788 --outfreq 10
os.system("python parcel.py --outfile "+str(outfile)+" --p_0 "+str(p_0)+" --T_0 "+str(T_0)+" --r_0 "+str(r_0)+" --w "+str(w)+" --dt "+str(dt)+" --sd_conc "+str(sd_conc)+" --z_max "+str(z_max)+" --z_min "+str(z_min)+" --top_stop "+str(top_stop)+" --nt "+str(int(nt))+" --aerosol \'"+str(aerosol)+"\'"+\
" --aerosol_sizes \'"+str(aerosol_sizes)+"\'"+\
"  --out_bin \'"+str(out_bin)+"\' --outfreq "+str(int(outfreq)))

