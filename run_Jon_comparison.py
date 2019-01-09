import os

outfile = 'plots/Jon_GCCN/Jon_GCCN.nc'
p_0     = 93844.6 
T_0     = 284.3 
#RH_0    = 0.85785275093866054 
r_0     = 0.007617065 
w       = 0.
dt      = .05
sd_conc = 10000 
z_max   = 0 #788

outfreq = 10
nt      = 10
outfreq = outfreq / dt
nt      = nt / dt

p_stp = 101325
T_stp = 273.25 + 15 

n_0 = 125e6
n_1 = 65e6

n_0_stp_corrected = n_0  * p_stp / p_0 * T_0 / T_stp
n_1_stp_corrected = n_1  * p_stp / p_0 * T_0 / T_stp

aerosol = '{"NaCl": {"kappa": 1.28, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": ['+str(n_0_stp_corrected)+','+str(n_1_stp_corrected)+']}}'

out_bin = '{"radii": {"rght": 1e-6, "moms": [0], "drwt": "dry", "nbin": 500, "lnli": "log", "left": 1e-9},'+\
          '"cloud": {"rght": 40e-6, "moms": [0], "drwt": "wet", "nbin": 500, "lnli": "log", "left": 1e-9}}'

#python parcel.py --outfile plots/Jon_GCCN/Jon_GCCN.nc --p_0 93844.6 --T_0 284.3 --RH_0 0.85785275093866054 --w 0. --dt .05 --sd_conc 1000 --z_max 788 --outfreq 10
os.system("python parcel.py --outfile "+str(outfile)+" --p_0 "+str(p_0)+" --T_0 "+str(T_0)+" --r_0 "+str(r_0)+" --w "+str(w)+" --dt "+str(dt)+" --sd_conc "+str(sd_conc)+" --z_max "+str(z_max)+" --nt "+str(int(nt))+" --aerosol \'"+str(aerosol)+"\' --out_bin \'"+str(out_bin)+"\'")
