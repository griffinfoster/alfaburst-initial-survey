#!/usr/bin/python

# calcalphamin.py
# Calculate the lower limit on FRB spectral index, given the standard candle
# model.
#
# Created by Jayanth Chennamangalam on 2014.10.25


"""Usage: calcalphamin.py [options]

Options:
    -h --help                         Display this usage information
    -f --f-low <freq>                 Lower limit of the frequency band in MHz
                                      [default: 142.84530816]
    -F --f-high <freq>                Upper limit of the frequency band in MHz
                                      [default: 148.9]
    -s --s-min <flux>                 Survey sensitivity in Jy
                                      [default: 61.75]
"""

import sys
#from docopt import docopt
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid.inset_locator as il


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


def calcD(z):
    c = 299792.458      # km s^-1
    H0 = 68.0           # km s^-1 Mpc^-1
    Om = 0.32
    Olambda = 0.68

    dz1 = 0.00001
    z1 = np.arange(0.0, z, dz1)
    # compute D in Mpc
    D = ((c * dz1)/H0) * np.sum(1.0/np.sqrt((Om * np.power(1 + z1, 3)) + Olambda))
    # convert to m
    D = D * 3.086e22

    return D


def calcL(alpha, fMin, fMax):
    SpeakTemp = 1e-19   # 1 Jy = 10^-26 W... = 10^-19 erg s^-1...
    zRef = 0.75
    fLoT13 = 1182e6     # 1182 MHz, in Hz
    fHiT13 = 1582e6     # 1582 MHz, in Hz

    L = ((fHiT13 - fLoT13) / (fHiT13**(alpha+1) - fLoT13**(alpha+1)))         \
        * SpeakTemp * 4 * np.pi * (calcD(zRef))**2                            \
        * (fMin**(alpha+1) - fMax**(alpha+1)) / ((1 + zRef)**(alpha-1))

    return L


def calcSpeak(L, z, alpha, fHi, fLo, fMin, fMax):
    Speak = (L * (1 + z)**(alpha-1)                                           \
             / (4 * np.pi * (calcD(z))**2                                     \
                * (fMin**(alpha+1) - fMax**(alpha+1))))                       \
            * ((fHi**(alpha+1) - fLo**(alpha+1)) / (fHi - fLo))
    # convert to Jy
    Speak /= 1e-19

    return Speak


def nearestIdx(array, value):                                                  
    idx = (np.abs(array - value)).argmin()                                      
    return idx


# input arguments
#args = docopt(__doc__, version="1.0.0")
#fLo = float(args["--f-low"])        # lower limit of the frequency band in MHz
#fHi = float(args["--f-high"])       # upper limit of the frequency band in MHz
#SMin = float(args["--s-min"])       # survey sensitivity limit in Jy

from optparse import OptionParser
o = OptionParser()
o.set_usage('%prog [options]')
o.set_description(__doc__)
o.add_option('--flow', dest='flow', type='float', default=142.84530816,
    help='lower limit of the frequency band in MHz')
o.add_option('--fhigh', dest='fhigh', type='float', default=148.9,
    help='upper limit of the frequency band in mhz')
o.add_option('--smin', dest='smin', type='float', default=61.75,
    help='survey sensitivity limit in Jy')
opts, args = o.parse_args(sys.argv[1:])

fLo = opts.flow
fHi = opts.fhigh
SMin = opts.smin

# convert to Hz
fLo = fLo * 1e6
fHi = fHi * 1e6

# constants
fMin = 1e10                 # 10 GHz, in Hz
fMax = 1e7                  # 10 MHz, in Hz
rate = 29

#zAvg = 0.177
zAvg = 0.132

#z = np.logspace(np.log10(0.0001), np.log10(0.11), 100)
#z = np.logspace(np.log10(0.0001), np.log10(zAvg), 100)
#z = np.logspace(np.log10(0.0001), np.log10(1.0), 100)
z = np.linspace(0.0001, 1.0, 1000)
#alpha = [1.2, 0.56, 0.44, 0.27, -0.1]
alpha = [0.19, 0.08, -0.11, -0.21]
alphaLabels = [0.2, 0.1, -0.1, -0.2]
pltStyle = ["g-.", "r-", "b:", "m--"]
Speak = np.zeros((len(alpha), len(z)))

# initialize tex-related stuff, and specify a font size
texInit(18)

plt.figure(2)#, figsize=(6, 8)) 
#plt.subplot(212)
for i in range(len(alpha)):
    L = calcL(alpha[i], fMin, fMax)
    # for each alpha, calculate Speak(z)
    for j in range(len(z)):
        Speak[i,j] = calcSpeak(L, z[j], alpha[i], fHi, fLo, fMin, fMax)
        #Speak[i,j] *= (2**0.4)
    label = r"$\alpha=%.1g$" % alphaLabels[i]
    plt.semilogy(z, Speak[i], pltStyle[i], label=label)

plt.axhline(y=SMin, xmin=0.0, xmax=1.0, color="k", linestyle="--")
plt.axvline(x=z[91], ymin=0.0, ymax=0.5, color="r", linestyle=":")
plt.axvline(x=z[78], ymin=0.0, ymax=0.5, color="g", linestyle=":")
plt.axvline(x=z[116], ymin=0.0, ymax=0.5, color="b", linestyle=":")
#plt.xlim(0.0, 0.11)
plt.xlim(0.0, zAvg)
#plt.xlim(0.0, 0.087)
plt.ylim(1e-2, 1e8)
#plt.annotate(r"${\rm (b)}$", xy=(0.03, 0.93), xycoords="axes fraction")
plt.xlabel(r"$z$")
plt.ylabel(r"$S_{\rm peak}~{\rm (Jy)}$")
plt.legend(loc="best", prop={"size":14})

plt.figure(1)
RTh = 10**4                             # 10^4 FRBs sky^-1 day^-1 above S_lim
RThLo = 0.5 * 10**4                     # 10^4 FRBs sky^-1 day^-1 above S_lim
RThHi = 1.6 * 10**4                     # 10^4 FRBs sky^-1 day^-1 above S_lim
#z = np.linspace(0.0, 1.0, 1000)
R = np.zeros(len(z))
RHi = np.zeros(len(z))
RLo = np.zeros(len(z))
for i in range(len(z)):
    D = calcD(z[i]) / 3.086e22          # Mpc
    Dlim = calcD(0.75) / 3.086e22       # Mpc
    R[i] = RTh * (D / Dlim)**3
    RHi[i] = RThHi * (D / Dlim)**3
    RLo[i] = RThLo * (D / Dlim)**3
plt.semilogy(z, R, "w")
plt.fill_between(z, RLo, RHi, alpha=1.0, facecolor="blue", edgecolor="white", hatch="\\\\\\")
plt.axhline(y=rate, xmin=0.0, xmax=0.25, color="k", linestyle=":")
plt.axvline(x=z[91], ymin=0.0, ymax=0.5, color="r", linestyle=":")
plt.axvline(x=z[78], ymin=0.0, ymax=0.5, color="g", linestyle=":")
plt.axvline(x=z[116], ymin=0.0, ymax=0.5, color="b", linestyle=":")
#plt.semilogy(z, RHi, "b--")
#plt.semilogy(z, RLo, "b--")
plt.xlim(0.0, 1.0)
plt.ylim(1.0, 1e5)
# plot ARTEMIS rate
# switched uplims and lolims to get the arrows right
plt.errorbar(zAvg, rate, yerr=10.0, capsize=3, fmt=None, ecolor="r",            \
             uplims=False, lolims=True)
# plot Thornton et al. (2013) rate
plt.plot(0.75, 1e4, "wx", markersize=10, markeredgewidth=4)
#plt.annotate(r"${\rm (a)}$", xy=(0.03, 0.93), xycoords="axes fraction")
plt.xlabel(r"$z$")
plt.ylabel(r"$R(< z)~({\rm bursts~sky}^{-1}{\rm ~day}^{-1})$")
print "nearestIdx = ", nearestIdx(R, rate)
print "z = ", z[nearestIdx(R, rate)]
print "Hi"
print "nearestIdx = ", nearestIdx(RHi, rate)
print "z = ", z[nearestIdx(RHi, rate)]
print "Lo"
print "nearestIdx = ", nearestIdx(RLo, rate)
print "z = ", z[nearestIdx(RLo, rate)]


'''
# show inset
ia = il.inset_axes(ax, width="40%", height="40%", loc=4)
ia.semilogy(z[:150], R[:150])
# plot ARTEMIS rate
# switched uplims and lolims to get the arrows right
plt.errorbar(0.11, rate, yerr=6.0, capsize=3, fmt=None, ecolor="r",            \
             uplims=False, lolims=True)
ia.set_xlim(0.05, 0.15)
ia.set_ylim(10.0, 100.0)
ia.set_xticks([0.05, 0.10])
ia.xaxis.tick_top()
'''

plt.tight_layout()
plt.show()

