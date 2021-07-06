import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from math import exp,sqrt,log
from scipy import optimize



A = 1.5e-9;
k = 0.61;
S = -0.02;


rw = np.logspace(-2,1,1000)*1e-6;
rd = pow((A/k*rw**2 -S/k*rw**3 ),1/3)

rd_prime = (2*A/3/k*rw - S/k*rw**2)*pow((A/k*rw**2 - S/k*rw**3),-2/3);
rd_prime_num = np.diff(rd)/np.diff(rw);

sigma = log(1.58);
mu = log(0.075);
Pd= np.zeros(len(rd))
Pd = 1/sigma/np.sqrt(2*np.pi)/(rd*1e6)*np.exp(-0.5*((np.log((rd*1e6))-mu)/sigma)**2);

# semilogx(rd*1e6,Pd)
Pw = Pd*rd_prime;
plt.plot(rw*1e6,Pw)
plt.plot(rd*1e6,Pd)
plt.xlabel('Radius [\mum]')
plt.ylabel('PDF')
plt.show()

# trapz(rd*1e6,Pd)
# trapz(rw*1e6,Pw)
