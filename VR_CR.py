#!/usr/bin/env python

from optparse import OptionParser
import sys
import os
import numpy as N
import datetime as DT 
import netCDF4 as ncdf
import glob
import pylab as P
from Plotting.cbook2 import nice_mxmnintvl, nice_clevels
from datetime import datetime

import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from scipy import ndimage

_bin_delta = 1
zbins = 2000. * N.arange(6)
radar_hgt = 373.

_prior_files = "Prior_*"
_plotfilename = "VrConsistencyRatio"

_cmin, _cmax = 0.0, 1.1   # limits for the single line plots...

variable = "VR"   # valid variables:  VR, DBZ
variable = 11   # valid variables:  VR, DBZ

sec_utime = "seconds since 1970-01-01 00:00:00"
auto_clevels = True

# definitions for the plot layout with multiple panels
left, width = 0.1, 0.5
bottom, height = 0.1, 0.5
bottom_h = left_h = left+width+0.03

rectC = [left, bottom, width, height]
rectX = [left, bottom_h, width, 0.2]
rectY = [left_h, bottom, 0.2, height]

kernel = N.array([[1,2,1],
                  [-2,0,-2],
                  [1,-2,1]])

def smfnc(x):
    return ndimage.gaussian_filter(x, sigma=0.75)
    
# ------------------------------------------------------------------------------------------
# Search functions  each returns a boolean array of T/F with the same dimensions as field

def getIndexVariable(field, variable):
    return (field == variable)
    
def getIndexGreaterThan(field, value):
    return ( field > value )
    
def getIndexLessThan(field, value):
    return ( field < value )
    
def getIndexGreaterThanOrEqual(field, value):
    return ( field >= value )
    
def getIndexLessThanOrEqual(field, value):
    return ( field <= value )

def getIndexEqual(field, value):
    return ( field == value )

def getIndexNotEqual(field, value):
    return ( field != value )

#####################################################################################################
# Main program

if __name__ == "__main__":

    print("\n<<<<<===========================================================================================>>>>>>")
    print("\n\n    CONSISTENCY RATIO PLOT       \n\n")

    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)

    parser.add_option("-d",  "--dir",    dest="dir",  type="string", help="Name of directory where Prior files are")
    parser.add_option("-t",  "--title",  dest="title", type="string", help="Name of plot and used as the name of the outputfile")
    parser.add_option(       "--noshow",  dest="noshow", default=False, action="store_true", help="Turn off screen plotting")

    (options, args) = parser.parse_args()

    if options.dir:
      dirname = os.path.join(options.dir,_prior_files)
    else:
      print("\n  ====>  No directory supplied, using %s as prefix \n" % _prior_files)
      dirname = os.path.join("./",_prior_files)

    if options.title:
      plotfilename = options.title
    else:
      plotfilename = _plotfilename

    file_list = glob.glob(dirname)
    file_list = sorted(file_list,key=os.path.getmtime)
    print(file_list)

    print("\nFirst file:  %s" % file_list[0])
    print("Last file:   %s\n" % file_list[-1])

    bin_delta = _bin_delta

    nbins = len(file_list) // bin_delta
    CR_TZ = N.zeros((zbins.size,nbins))
    CR_T  = N.zeros((nbins))
    CR_Z  = N.zeros((zbins.size))
    m = -1
    
    fig = P.figure(figsize=(12,12))
    fig.text(0.72, 0.72, "\n\nConsistency\n\nRatio", size=20, va="baseline", ha="center", multialignment="center")
    fig.text(0.72, 0.65, plotfilename, size=16, va="baseline", ha="center", multialignment="center")

    print("\n 2D CONSISTENCY RATIO CALCULATIONS......\n")
    
    datebins = []
    secsbins = []

    for n, file in enumerate(file_list):

        f = ncdf.Dataset(file)

        if n % bin_delta == 0:
        
            HxfL     = []
            depL     = []
            errorL   = []
            secsL    = []
            zL       = []
            kindL    = []

        HxfL.append(f.variables["Hxf"][:])
        zL.append(f.variables["z"][:] - radar_hgt)
        depL.append( f.variables["value"][:] - f.variables["Hxfbar"][:] )
        secsL.append(f.variables["secs"][:])
        errorL.append(N.sqrt(f.variables["error"][:]))
#       kindL.append(f.variables["type"][:])
        kindL.append(f.variables["kind"][:])

        if n % bin_delta == bin_delta-1:
            Hxf   = N.concatenate(HxfL, axis=0)
            z     = N.concatenate(zL, axis=0)
            dep   = N.concatenate(depL, axis=0)
            secs  = N.concatenate(secsL, axis=0)
            error = N.concatenate(errorL, axis=0)
            kind  = N.concatenate(kindL, axis=0)
            datebins.append((ncdf.num2date(secs[1],units=sec_utime)).strftime("%Y%m%d%H%M%S"))

            m = m + 1

            index_kind  = getIndexVariable(kind,variable)
       
            for k in N.arange(zbins.size-1):
                index1                = getIndexGreaterThanOrEqual(z, zbins[k])
                index2                = getIndexLessThan(z, zbins[k+1])
                index                 = index_kind & index1 & index2
                if N.sum(index == True) > 2: 
                    d           = dep[index]
                    obs_var     = error[index]
                    Hxftmp      = Hxf[index,:]
                    Hxf_var     = Hxftmp.var(ddof=1, axis=1).mean()
                    inno_var    = N.mean((d - d.mean())**2)
                    consi_ratio = (obs_var[1]**2 + Hxf_var) / inno_var
                    print("%s  NOBS: %5.5d    %3.3s: %3.1f  ZBIN:  %f  %f  RMSI: %6.3f  M-Innov: %7.3f  Spread: %6.3f  CRatio: %7.4f " \
                    % (file[-22:-3], d.size, "VR", obs_var[1], zbins[k], zbins[k+1], \
                    N.sqrt(inno_var), d.mean(), N.sqrt(obs_var[1]**2 + Hxf_var), consi_ratio))
                    CR_TZ[k,m] = consi_ratio
          
        f.close()

# Plotting

    axC = P.axes(rectC)

    zmin, zmax, cint = nice_mxmnintvl(zbins.min()/1000., zbins.max()/1000., outside=True, cint=1.0)
    zmin = 0.0

    if auto_clevels:
        cmin, cmax, cint, clevels = nice_clevels(0, 3, outside=False, cint = 0.1)
    else:
        clevels = N.arange(spread_limits[0], spread_limits[1], spread_limits[2])

    cs1=axC.contourf(datebins, zbins/1000., CR_TZ, clevels, cmap=cm.get_cmap('YlOrRd'))
    cs2=axC.contour(datebins,  zbins/1000., CR_TZ, cs1.levels, colors='k')

    start = datebins[0]
    end   = datebins[-1]
    s     = datetime.strptime(start, "%Y%m%d%H%M%S")
    e     = datetime.strptime(end, "%Y%m%d%H%M%S")

    axC.set_xlim(start, end)
    axC.set_ylim(zmin,zmax)

    maj_loc = mdates.MinuteLocator(interval=2)
    axC.xaxis.set_major_locator(maj_loc)
    dateFmt = mdates.DateFormatter('%H:%M')
    axC.xaxis.set_major_formatter(dateFmt)

    min_loc   = mdates.MinuteLocator(interval=2)
    axC.xaxis.set_minor_locator(min_loc)

    labels = axC.get_xticklabels()
    P.setp(labels, rotation=40, fontsize=10)

    axC.clabel(cs2, inline=1, fontsize=10, fmt="%1.1f")
    axC.set_ylabel("Height (km)")
    axC.set_xlabel("Time")

#===================================================================================================
# time series Consistency Ratio

    print("\n TIME-SERIES CONSISTENCY RATIO CALCULATIONS......\n")

    datebins = []
    m = -1

    for n, file in enumerate(file_list):

        f = ncdf.Dataset(file)

        if n % bin_delta == 0:
        
            HxfL     = []
            depL     = []
            errorL   = []
            secsL    = []
            zL       = []
            kindL    = []

        HxfL.append(f.variables["Hxf"][:])
        zL.append(f.variables["z"][:] - radar_hgt)
        depL.append( f.variables["value"][:] - f.variables["Hxfbar"][:] )
        secsL.append(f.variables["secs"][:])
        errorL.append(N.sqrt(f.variables["error"][:]))
#       kindL.append(f.variables["type"][:])
        kindL.append(f.variables["kind"][:])

        if n % bin_delta == bin_delta-1:
            Hxf   = N.concatenate(HxfL, axis=0)
            z     = N.concatenate(zL, axis=0)
            dep   = N.concatenate(depL, axis=0)
            secs  = N.concatenate(secsL, axis=0)
            error = N.concatenate(errorL, axis=0)
            kind  = N.concatenate(kindL, axis=0)
            datebins.append((ncdf.num2date(secs[1],units=sec_utime)).strftime("%Y%m%d%H%M%S"))

            m = m + 1

            index = getIndexVariable(kind,variable)
       
            if N.sum(index == True) > 2: 
                d           = dep[index]
                obs_var     = error[index]
                Hxftmp      = Hxf[index,:]
                Hxf_var     = Hxftmp.var(ddof=1, axis=1).mean()
                inno_var    = N.mean((d - d.mean())**2)
                consi_ratio = (obs_var[1]**2 + Hxf_var) / inno_var
                print("%s  NOBS: %5.5d    %3.3s: %3.1f  ZBIN:  %f  %f  RMSI: %6.3f  M-Innov: %7.3f  Spread: %6.3f  CRatio: %7.4f " \
                % (file[-22:-3], d.size, "VR", obs_var[1], 0.0, zbins.max(), \
                N.sqrt(inno_var), d.mean(), N.sqrt(obs_var[1]**2 + Hxf_var), consi_ratio))
                CR_T[m] = consi_ratio
      
        f.close()

# Plotting

    axX = P.axes(rectX)
#   cmin, cmax, cint = nice_mxmnintvl(CR_T.min(), CR_T.max(), outside=True, cint=1.0)  # Use this to get limits of the plot
    s = datebins[0]
    e = datebins[-1]
    axX.plot(datebins, CR_T, lw=2.0, color='k')
    axX.set_xlim(s, e)
    axX.set_ylim(_cmin, _cmax)
    axX.set_xticklabels([])
    axX.set_ylabel("Consistency Ratio")
    axX.grid(True)

    maj_loc = mdates.MinuteLocator(interval=10)

#===================================================================================================
# Height Consistency Ratio
    
    print("\n ZBIN-D CONSISTENCY RATIO CALCULATIONS......\n")

    datebins = []
    m = -1
        
    HxfL     = []
    depL     = []
    errorL   = []
    secsL    = []
    zL       = []
    kindL    = []

    for n, file in enumerate(file_list):

        f = ncdf.Dataset(file)

        HxfL.append(f.variables["Hxf"][:])
        zL.append(f.variables["z"][:] - radar_hgt)
        depL.append( f.variables["value"][:] - f.variables["Hxfbar"][:] )
        secsL.append(f.variables["secs"][:])
        errorL.append(N.sqrt(f.variables["error"][:]))
#       kindL.append(f.variables["type"][:])
        kindL.append(f.variables["kind"][:])

        f.close()

    Hxf   = N.concatenate(HxfL, axis=0)
    z     = N.concatenate(zL, axis=0)
    dep   = N.concatenate(depL, axis=0)
    secs  = N.concatenate(secsL, axis=0)
    error = N.concatenate(errorL, axis=0)
    kind  = N.concatenate(kindL, axis=0)

    index_kind  = getIndexVariable(kind,variable)

    for k in N.arange(zbins.size-1):
        index1                = getIndexGreaterThanOrEqual(z, zbins[k])
        index2                = getIndexLessThan(z, zbins[k+1])
        index                 = index_kind & index1 & index2
        if N.sum(index == True) > 2: 
            d           = dep[index]
            obs_var     = error[index]
            Hxftmp      = Hxf[index,:]
            Hxf_var     = Hxftmp.var(ddof=1, axis=1).mean()
            inno_var    = N.mean((d - d.mean())**2)
            consi_ratio = (obs_var[1]**2 + Hxf_var) / inno_var
            print("%s  NOBS: %5.5d    %3.3s: %3.1f  ZBIN:  %05.0f  %05.5f  RMSI: %6.3f  M-Innov: %7.3f  Spread: %6.3f  CRatio: %7.4f " \
            % (file[-22:-3], d.size, "VR", obs_var[1], zbins[k], zbins[k+1], \
            N.sqrt(inno_var), d.mean(), N.sqrt(obs_var[1]**2 + Hxf_var), consi_ratio))
            CR_Z[k] = consi_ratio
    
# Plotting
  
    axY = P.axes(rectY)
#   cmin, cmax, cint = nice_mxmnintvl(CR_Z.min(), CR_Z.max(), outside=True, cint=1.0)  # Use this to get limits of the plot
    
    axY.plot(CR_Z, zbins/1000., lw=2.0, color='k')
    axY.set_ylim(0.0,zbins.max()/1000.)
    axY.set_xlim(_cmin,_cmax)
    axY.set_yticklabels([])
    axY.set_xlabel("Consistency Ratio")
    axY.grid(True)
    
#===================================================================================================
# SAVE FILE AND SHOW PLOT

    P.savefig(plotfilename+".pdf")
    if not options.noshow:
        P.show()

