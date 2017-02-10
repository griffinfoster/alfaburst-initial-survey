#!/usr/bin/python

# optdm.py
# Compute optimal DM step size
#
# Created by Jayanth Chennamangalam
# 2014.10.19

import numpy as np
import scipy.special as ss
import matplotlib as mp
import matplotlib.pyplot as plt


# function definitions
def texInit(fontsize):
    # set plotting font properties
    font = {"family" : "serif",
            "weight" : "regular",
            "size"   : fontsize}
    plt.rc("font", **font)
    # force matplotlib to use Type 1 fonts instead of Type 3
    mp.rcParams["ps.useafm"] = True
    mp.rcParams["pdf.use14corefonts"] = True
    mp.rcParams["text.usetex"] = True


texInit(16)

#f = 145                 # MHz
#deltaf = 6              # MHz
f = 1375                 # MHz
deltaf = 56              # MHz
a = -6.46               # Bhat et al. 2004
#a = -9.5                # Lorimer et al. 2007
b = 0.154
c = 1.07
alpha = 3.86
ts = 0.32768            # ms
#Wint = np.arange(ts, 64 * ts, ts)
Wint = [2, 5]    # ms
pltStyle = ["b--", "r-"]

#DM = np.linspace(1, 320, 219)
DM = np.linspace(1, 1500, 500)
tau = 10**(a
           + (b * np.log10(DM))
           + (c * (np.log10(DM)**2))
           - (alpha * np.log10(f * 1e-3)))
tau = 0.1 * tau

#deltaDM = 0.1
deltaDM = 1.

fGHz = f * 1e-3

for i in range(len(Wint)):
    zeta = 6.91e-3 * deltaDM * deltaf / (np.sqrt(Wint[i]**2 + tau**2) * fGHz**3) 
    S = np.sqrt(np.pi) * ss.erf(zeta) / (2 * zeta)
    plt.plot(DM, S, pltStyle[i], label=r"$%1.0f~{\rm ms}$" % Wint[i])

#plt.xlim(0, 320)
#plt.ylim(5*1e-2, 1e2)
plt.legend(loc="best")
plt.xlabel(r"${\rm DM~(cm}^{-3}{\rm ~pc)}$")
plt.ylabel(r"$S(\delta{\rm DM})/S$")
plt.show()

