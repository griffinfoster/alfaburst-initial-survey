#!/usr/bin/python

#
# calcfrbvol.py
# Given a DM limit of an FRB search and a set of beam pointing in Az-El,
# compute the search volume of the LOFAR-Chilbolton survey
# Dependency: NE2001
#
# Jayanth Chennamangalam
# Parts based on code written by Evan Keane
# 2014.06.25
#

"""Usage: calcfrbvol.py [options]

Options:
    -h --help                         Display this usage information
    -d --dm <dm>                      Dispersion measure limit of search
                                      in cm^-3 pc
                                      [default: 320]
    -b --beams <filename>             File containing beam Az-El
                                      [default: "azel"]
    -g --graphics                     Turn on graphics

"""

import sys
#from docopt import docopt
import numpy as np
import ephem
import subprocess as sp
import matplotlib as mp
import matplotlib.pyplot as plt


def Dist(z):
    c = 299792.458      # km s^-1
    H0 = 68.0           # km s^-1 Mpc^-1
    Om = 0.32
    Olambda = 0.68

    dz1 = 0.001
    z1 = np.arange(0.0, z, dz1)
    # compute D in Mpc
    D = ((c * dz1)/H0) * np.sum(1.0/np.sqrt((Om * np.power(1 + z1, 3)) + Olambda))
    ## convert to m
    #D = D * 3.086e22
    # convert to Gpc
    D = D * 1e-3

    return D


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


# initialize tex-related stuff, and specify a font size                         
texInit(16)

#args = docopt(__doc__, version="1.0.0")
#
## input arguments
#DMLim = float(args["--dm"])         # DM limit in cm^-3 pc
#FileBeams = args["--beams"]         # File containing beam Az-El
#Plot = args["--graphics"]           # plot flag

from optparse import OptionParser
o = OptionParser()
o.set_usage('%prog [options]')
o.set_description(__doc__)
o.add_option('--dm', dest='dm', type='float', default=320.,
    help='DM limit in cm^-3 pc')
o.add_option('--beams', dest='beams', default=None,
    help='File containing beam Az-El')
o.add_option('--graphics', dest='graphics', action='store_true',
    help='plot flag')
opts, args = o.parse_args(sys.argv[1:])

DMLim = opts.dm
FileBeams = opts.beams
Plot = opts.graphics

# ASSUMPTIONS:
# 1. All FRBs originate outside the Galaxy
# 2. Some electron density in the IGM (Ioka 2003 and Inoue 2004)
# 3. Cosmological assumptions involved in d = (c * z) / H0 (?)

#Program = "/home/jayanth/programs/NE2001/bin.NE2001/NE2001"
Program = "/home/griffin/Documents/Manuscripts/alfaburst-suvery/artemisScripts/NE2001/bin.NE2001/NE2001" # HARDCODE
TailProcs = "| tr -s ' ' | sed 's/^[ \t]*//'"
# choose some large distance outside the Galaxy, so that NE2001 returns the
# DM for the entire line-of-sight within the Galaxy
DistBound = 100                 # kpc
# direction: distance to DM
DirDist2DM = "-1"
# direction: DM to distance
DirDM2Dist = "1"

# read beam info file and get azimuths and elevations of the LOFAR beams
AzEl = np.loadtxt("azel", dtype=float, delimiter=" ")
#AzEl = np.loadtxt("azeltest", dtype=float, delimiter=" ")

DMHost = 0.0

# constants
c = 299792458           # m s^-1
H0 = 68                 # km s^-1 Mpc^-1
UnitFactor = 1e-6       # to get the distance in Gpc

# observation start time
JDObs = 2455822.0               # JD
ObsDur = 1.0                    # duration of observation in days

# site details (HBA)
# taken from http://www.lofar-uk.org/wiki/uploads/Main/LOFAR-Chilbolton.pdf
Site = ephem.Observer()         # creates empty observer structure
Site.lat = str(51.143530)       # latitude, in degrees
Site.long = str(-1.434450)      # longitude, in degrees
Site.elevation = 0.0            # altitude in metres
JDRef = ephem.julian_date(0)    # reference Julian Day, 1899.12.31 noon
Site.date = JDObs - JDRef       # observation date
BeamSize = 1.93772018371        # measured by 0950, in degrees

# time taken for the sky to drift over one complete beam
# = ((1 day / 360 degrees) * BeamSize degrees) days
#TimeBeamScan = ((1.0 / 360) * BeamSize)
#print "Time for one beam scan = %g minutes" % (TimeBeamScan * 24 * 60)
# measured using 0950
TimeBeamScan = 7.72971747174    # minutes
print "Time for one beam scan = %g minutes" % TimeBeamScan
# convert to days
TimeBeamScan /= (24 * 60)

deg2rad = np.pi / 180
rad2deg = 180.0 / np.pi

if Plot:
    # set plotting font properties
    font = {"family" : "serif",
            "weight" : "regular",
            "size"   : 18}
    plt.rc("font", **font)
    # force matplotlib to use Type 1 fonts instead of Type 3
    mp.rcParams["ps.useafm"] = True
    mp.rcParams["pdf.use14corefonts"] = True
    mp.rcParams["text.usetex"] = True

    plt.subplot(111, projection="hammer")
    #####
    #DMMap = np.loadtxt("dm320.gal.out")
    #l = np.linspace(-np.pi, np.pi, 360)
    #b = np.linspace(-(np.pi / 2), (np.pi / 2), 181)
    #DMMap = DMMap.T
    #DMMap = np.array([np.roll(DMMap[i], DMMap.shape[1] / 2)                   \
    #                  for i in range(DMMap.shape[0])])
    #if l > 180.0:                                                                   
    #    lPlot = 360 - l                                                                 
    #else:                                                                           
    #    lPlot = -l  
    #plt.pcolormesh(lPlot, b, DMMap, cmap="RdYlBu")
    ##cbar = plt.colorbar(img, orientation="vertical")
    ##cbar.set_label("DM (cm^-3 pc)")
    ####

#GalLat = []
#GalLong = []
#RedShift = []

# number density of galaxies
n = 1e-2        # Mpc^-3
# number of galaxies
N = 0.0

# volume of universe
V = 0.0
# area of sky
A = 0.0
# average z
zAvg = 0.0
# maximum z
zMax = -1000.0
BeamID = 1
TotalPointings = 0.0
EGPointings = 0.0
for Az, El in AzEl:
    sys.stdout.write("Processing beam %d..." % BeamID)
    sys.stdout.flush()
    JDNow = JDObs
    tObs = 0.0
    while tObs < ObsDur:
        # increment elapsed observation time
        tObs += TimeBeamScan
        # update current time
        JDNow += TimeBeamScan
        Site.date = JDNow - JDRef
        # get RA and dec.
        RA, Dec = Site.radec_of(Az, El)
        # convert (RA, Dec) to Galactic co-ordinates
        GalPos = ephem.Galactic(ephem.Equatorial(RA, Dec))
        # the internal representation (see 'repr(GalPos.lat)') is in radians,
        # so convert to degrees
        b = GalPos.lat / deg2rad
        l = GalPos.lon / deg2rad

        # compute non-Galactic DM contribution
        # build NE2001 command string
        Cmd = " ".join((Program, str(l), str(b), str(DistBound), DirDist2DM,
                        TailProcs))
        # run NE2001 command
        Buf = sp.check_output(Cmd, shell=True, universal_newlines=True)
        # get the 7th line
        Cols = Buf.split("\n")[7].split(" ")
        DMGal = float(Cols[0])
        #DMGal = 0.0
        #DMHost = DMGal
        DMHost = 100.0
        DMDiff = DMLim - DMGal - DMHost
        TotalPointings += 1
        if DMDiff > 0:      # if outside the Galaxy, and able to be seen
                            # outside the host galaxy with our DMLim
            # using the relation DM ~ 1200 z cm^-3 pc for z <= 2 (Lorimer
            # et al. 2007), based on Ioka (2003) and Inoue (2004)
            EGPointings += 1
            z = DMDiff / 1200
            zAvg += z
            if z > zMax:
                zMax = z
            # compute the Hubble distance
            #D = (c * z * UnitFactor) / H0
            D = Dist(z)
            # volume, assuming a cone, V = pi r^2 h / 3, where r = d theta / 2,
            # and h = d
            # NOTE: this underestimates the volume as it does not consider the
            # gaps between the circular beams
            #V += np.pi * (BeamSize * deg2rad)**2 * D**3 / 12
            # volume, assuming a pyramid, V = l w h / 3, where l = w = d theta,
            # and h = d
            # NOTE: this overestimates the volume by a tiny amount (in the
            # first and last pointings of the scan)
            Vp = ((BeamSize * deg2rad)**2 * D**3) / 3
            Np = Vp * 1e9 * n
            N += Np
            #print Np
            V += ((BeamSize * deg2rad)**2 * D**3) / 3
            # area, treat as square, not circle, as this is a drift scan
            # (overestimates a bit, in the first and last pointings; only
            # relevant if DMHost = 0.0)
            A += (np.pi * (BeamSize / 2)**2)

        if Plot:
            # build NE2001 command string
            Cmd = " ".join((Program, str(l), str(b), str(DMLim), DirDM2Dist,
                            TailProcs))
            # run NE2001 command
            Buf = sp.check_output(Cmd, shell=True, universal_newlines=True)
            # get the 6th line
            Cols = Buf.split("\n")[6].split(" ")
            if l > 180.0:                                                                   
                lPlot = 360 - l                                                                 
            else:                                                                           
                lPlot = -l  
            if (">" == Cols[0]):
                plt.plot(lPlot * deg2rad, b * deg2rad, "go", markersize=5)
            else:
                plt.plot(lPlot * deg2rad, b * deg2rad, "ro", markersize=5)
    print "DONE"
    BeamID += 1

print "Total volume surveyed = %2.3e Gpc^3" % V, "= %2.3e Mpc^3" % (V * 1e9)
print "Total area surveyed = %2.3e sq. deg." % A
print "Maximum z = %2.3f" % zMax
print "Average z = %2.3f" % (zAvg / EGPointings)
print "Total number of galaxies (assuming n = %2.3f" % n, "Mpc^-3) = %2.3f" % N

print "Fraction of extragalactic pointings =", (EGPointings / TotalPointings)

if Plot:
    ticks, labels = plt.xticks()
    newTicks = ticks[::-1] * rad2deg
    for i in range(len(newTicks)):
        if newTicks[i] < 0.0:
            newTicks[i] = 360 + newTicks[i]
    plt.xticks(ticks, map(lambda val: r"$%2.0f$" % val, newTicks))
    plt.grid(True)
    plt.xlabel(r"$l$")
    plt.ylabel(r"$b$")
    plt.show()

