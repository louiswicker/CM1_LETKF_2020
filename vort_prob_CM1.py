#!/usr/bin/env python
#

# This python script plots vorticity swaths overlaid from different ensemble members
# produced from a pyencommas run.  Its intent is for an ensemble tornado prediction study.

# Import needed libraries

from optparse import OptionParser
import glob
import pickle
import time
import netCDF4 as ncdf
import os, sys
from numpy import *
import numpy as N
import datetime as DT
import matplotlib.pyplot as P
import matplotlib
import matplotlib.cm as cm
from matplotlib import ticker
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid import AxesGrid
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import ctables
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
from scipy import signal

from matplotlib import mpl
import cbook2 as cbook
import ctables
import ens
from matplotlib.colors import BoundaryNorm

from state_vector import vswath
from Plotting import shapefiles

import fsrc.cressman as cress

#===================================================================================================
def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = N.mgrid[-size:size+1, -sizey:sizey+1]
    g = N.exp(-(x**2/float(size) + y**2/float(sizey)))
    return g / g.sum()

#===================================================================================================
def gauss_prob(im, n, ny=None):
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """

    g = gauss_kern(n/2)

    nx     = im.shape[0]
    ny     = im.shape[1]
    buf    = n / 2
    improc = N.zeros(im.shape)

    improc[buf:nx-buf,buf:ny-buf] = signal.convolve2d(im, g, mode='valid')

    return(improc)

#===================================================================================================
def nearest_N_prob(im, n):
    """
       Uses a 2D convolution operator to see if any N x N neighbors have a non-zero entry
    """

    g = N.ones([n,n])
    g /= g.sum()

    nx = im.shape[0]
    ny = im.shape[1]
    buf = n / 2
    improc = N.zeros(im.shape)

    improc[buf:nx-buf,buf:ny-buf] = signal.convolve2d(im, g, mode='valid')

    return( N.where(improc > 0.0, 1.0, 0.0) )
    
#---------------------------------------------------------------------------------------------------
# PARAMETERS

_yes_plot     = True
_write_prob   = False  # set this to True if you want to write prob data out to a file

_plotpcolormesh = False
_plotfilename   = "par_2km"
_runname        = 'par_2km'                 # Runname prefix
_threshold      = 0.005                    # threshold of vorticity
_threshold2     = -0.010
member_start    = 1
member_end      = 40
member_step     = 1
ens_size        = (member_end - member_start + 1) / member_step        # Number of ensemble members
level           = 8                                                    # Vertical model level to plot
_neighbor       = 4

# _default_color_table = cm.jet
# _default_color_table = ctables.__getattribute__("Sobash")
_default_color_table = cm.OrRd
# _default_color_table = ctables.__getattribute__("Thomas")
# _default_color_table = ctables.__getattribute__("WoF")member_start   = 1

# Set domain boundaries
#-----------------------
# Trim out the plot...zoomed in view
xplot         = [-100000.,60000,4000.]
yplot         = [-50000.,110000.,4000.]
# Trim out the plot...zoomed in view
xplot         = [-120000.,20000,4000.]
yplot         = [-20000.,120000.,4000.]


# Stuff for overlays from DBZ files
#-----------------------
dbz_threshold = 40
_overlayfile = ["PAR_1min_DBZ_F2020_Valid_T+30min.nc"]
_overlayfile = ["PAR_5min_DBZ_F2020_Valid_T+10min.nc"]
_overlayfile = ["PAR_5min_DBZ_F2020_Valid_T+20min.nc"]
_overlayfile = ["PAR_5min_DBZ_F2020_Valid_T+30min.nc"]
_overlayfile = None

#===================================================================================================
def mtokm(val,pos):
  """Convert m to km for formatting axes tick labels"""
  val=val/1000.0
  return '%i' % val
  
#===================================================================================================
# Read U & V and compute return vorticity

def CM1_get_Wz(files, exper):

  xoffset   = exper['xoffset']
  yoffset   = exper['yoffset']

  wz_ens = N.zeros((ens_size,level,exper['cm1namelist'][0][2],exper['cm1namelist'][1][2]))
  
  for n, file in enumerate(files[0:ens_size]):
    f         = ncdf.Dataset(file, "r")
    if n == 0:
      xc      = f.variables['xh'][:] + xoffset
      yc      = f.variables['yh'][:] + yoffset
      
    wz_ens[n] = ens.ComputeWZ(xc, yc, f.variables['ua'][0:level], f.variables['va'][0:level])
  
  return wz_ens, xc, yc

#===================================================================================================
# Main function defined to return correct sys.exit() calls

fig_type = 'pdf'                         # Figure type (pdf, png, eps, etc.)

print
print "<<<<<===========================================================================================>>>>>>"
print
print "                                        "
print
print "                 VORT SWATH PLOT         "
print

usage = "usage: %prog [options] arg"
parser = OptionParser(usage)

parser.add_option("-d",  "--dir",    dest="dir",      default=None, type="string", help = "Experiment directory - looks for a *.exp file in this directory")
parser.add_option("-e",  "--exp",    dest="exp",      default=None, type="string", help = "Path to the pickled experiment file generated by run")
parser.add_option("-t",  "--time",   dest="datetime", default=None, type= "string", help = "Usage:  --time 2003,5,8,21,0,0")
parser.add_option(       "--fcst",   dest="fcst",     default=None, type="int", nargs=2, help="Length of time, in model secs for swatch and number of sec between analysis times")
parser.add_option(       "--title",  dest="title",    default=None, type="string", help="name of plot")
parser.add_option(       "--noshow", dest="noshow",   default=False, action="store_true", help="Turn off screen plotting")

(options, args) = parser.parse_args()

if options.dir and options.exp == None:
    options.exp = glob.glob(os.path.join(options.dir, "*.exp" ))[0]
    print("\n ==> ENS_MAIN: found experiment files %s" % options.exp)

if options.exp == None:
    parser.print_help()
    print("\n ==> ENS_MAIN: ERROR --> Experiment's filename not supplied..EXITING!!!")
    sys.exit(-1)
else:
    with open(options.exp, 'rb') as f:
        exper = pickle.load(f)
    if options.dir:
      exper['base_dir'] = options.dir

if options.datetime == None:
    parser.print_help()
    print("\n ==> ENS_MAIN: ERROR --> time for files not supplied..EXITING!!!")
    sys.exit(-1)
else:
    dt = [int(i) for i in options.datetime.split(',')]
    myDT = DT.datetime(dt[0],dt[1],dt[2],dt[3],dt[4],dt[5])

if options.title:
    plotfilename = options.title
    if plotfilename.find("VP") == -1:
      plotfilename = plotfilename+"_VP"
else:
    plotfilename = options.dir+"_VP"

if options.fcst == None:
    parser.print_help()
    print("\n ==> ENS_MAIN: ERROR -->forecast time for files and secs between forecast times not supplied..EXITING!!!")
    sys.exit(-1)
else:
     fcstlen, dt_files = options.fcst[0], options.fcst[1]
     ntimes =  1 + fcstlen/dt_files

print("\nStarting time of the forecast:     %s" % myDT.strftime("%H:%M:%S"))
print("\nEndding time of the forecast:      %d" % fcstlen)
print("\nTime between files to be read in:  %d" % dt_files)
print("\nNumber of files to be read in:     %d" % ntimes)
print("\nSize of the ensemble be read in:   %d" % ens_size)

# Create masterlist of files and date and times

master_list = []
master_DT   = []

for n in N.arange(ntimes):
    master_DT.append(myDT + DT.timedelta(seconds=int(n*dt_files)))
    master_list.append(ens.FindRestartFiles(exper, master_DT[n], ret_exp=False, ret_DT=False))

ny, nx = exper['cm1namelist'][1][2],exper['cm1namelist'][0][2]
dx, dy = exper['cm1namelist'][3][2],exper['cm1namelist'][4][2]
lat0, lon0 = exper['lat0'], exper['lon0']

vvort  = N.zeros((ntimes,ens_size,level,ny,nx))
vvmean = N.zeros((ntimes,ens_size,ny,nx))

for n in N.arange(ntimes):
  print("Reading and computing Vert. Vort for ensemble at time %s" % master_DT[n].strftime("%H:%M:%S"))
  vvort[n], xc, yc = CM1_get_Wz(master_list[n], exper)
  vvmean[n]        = vvort[n].mean(axis=1)
  
# Coming out of this, you now have VVort(ntimes, nens, 0:level, ny, nx), VVmean[ntimes,nens, ny, nx], xx[ntimes,ny,nx], yy[ntimes,ny,nx]

if xplot:
  x_swath = xplot[0] + xplot[2]*N.arange(1+N.int((xplot[1]-xplot[0])/xplot[2]))
  y_swath = yplot[0] + yplot[2]*N.arange(1+N.int((yplot[1]-yplot[0])/yplot[2]))
else:
  x_swath = xc
  y_swath = yc
  
vort   = N.zeros((x_swath.size,y_swath.size), order='F')
xx, yy = N.meshgrid(xc, yc)
  
for n in N.arange(ens_size):

  for m in N.arange(ntimes):
  
    if m == 0:
      xobs   = xx.flatten()
      yobs   = yy.flatten()
      obs    = vvmean[m,n].flatten()
      xobs2  = xx.flatten()
      yobs2  = yy.flatten()
      obs2   = vvmean[m,n].flatten()
    else:
      xobs   = N.concatenate((xobs,xx.flatten()))
      yobs   = N.concatenate((yobs,yy.flatten()))
      obs    = N.concatenate((obs,vvmean[m,n].flatten()))
      xobs2  = N.concatenate((xobs,xx.flatten()))
      yobs2  = N.concatenate((yobs,yy.flatten()))
      obs2   = N.concatenate((obs,vvmean[m,n].flatten()))

  mask1 = where(obs > _threshold)
  obs[mask1] = 1.0

  if obs[mask1].size > 0:
    vort = vort + cress.cressman(xobs[mask1], yobs[mask1], obs[mask1], x_swath, y_swath, _neighbor*dx)

prob1 = (vort.transpose()/ens_size)*100.
# 
print "\nMax probability is:  ",prob1.max(),"\n"
print "\nMin probability is:  ",prob1.min(),"\n"

# if obs2.size > 0:
#   print "\nNumber of obs > %4.4f is %d\n" % (_threshold2, obs2.size)
  

#-----------------
if _yes_plot:

  fig, ax = P.subplots(1, 1, figsize=(12,10), sharex=True, sharey=True)

# Create map coordinates
#-------------------

  map = ens.mymap(x_swath, y_swath, lat0, lon0, ax = ax, shape_env=shapefiles, counties=True, noticks=False)
  lon2d, lat2d, x2d, y2d = map.makegrid(x_swath.size, y_swath.size, returnxy=True)

  clevels  = [0.0, 9., 10., 20, 30., 40., 50., 60., 70., 80., 90., 100.]
  clevels  = [0.0, 25., 30., 40.0, 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.]
  clevels  = [1., 10., 20., 40., 60., 80., 100.]

  cmap     = _default_color_table
  plot     = map.contourf(x2d, y2d, prob1, cmap=cmap, levels=clevels, ax=ax)
  cbar     = map.colorbar(plot,location='bottom',pad="5%")
  cbar.set_label('Ensemble Mean Probability of Vert. Vort > %4.4f %s' % (_threshold, "$s^{-1}$"))
    
# plot     = map.contour(xx, yy, prob1.transpose(), colors='k', alpha=0.2, levels=[10.])
# plot     = map.contour(xx, yy, prob1.transpose(), colors='k', alpha=0.3, levels=[50.])
# plot     = map.contour(xx, yy, prob1.transpose(), colors='k', alpha=1.0, levels=[90.])

  if _threshold2 > 0.0:
    vlat, vlon = cbook.dxy_2_dll(xobs2, yobs2, lat0, lon0, degrees=True)
    xx2, yy2 = map(vlon, vlat)
    map.plot(xx2, yy2, 'o', markersize=2, markerfacecolor='blue')
    

# OVERLAY dBZ information
#-------------------------
  if _overlayfile != None:

    fover = netcdf.Dataset(overlayfile, "r")
    prob_dbz = fover.variables['PROB'][...]
    mean_dbz = fover.variables['MEAN'][...]
    obs_dbz  = fover.variables['OBS'][...]

  # Plot mean threshold contour (using negatives to automatically get dashed contour)
    plot     = map.contour(xx, yy, mean_dbz, colors='r', linewidths=2.0, levels=[-dbz_threshold])
  
  # Plot Actual thresholded reflectivity from radar
    plot     = map.contour(xx, yy, obs_dbz, colors='r', linewidths=4.0, levels=[dbz_threshold])
    plot     = map.contourf(xx, yy, obs_dbz, colors='r', alpha = 0.1, levels=[dbz_threshold,80.])
    fover.close()

# Label some things
#------------------

  DTstart      = master_DT[0]
  DTend        = master_DT[-1]
  DT_yymmdd    = DTstart.strftime("%m-%d-%Y")
  DT_hhmmss    = DTstart.strftime("%H:%M")
  DTend_hhmmss = DTend.strftime("%H:%M")

  title = "%s\n%d minute FCST (%s to %s UTC)\n Mean Vertical Vorticity computed in layer below %3.2f KM" \
      % (DT_yymmdd, fcstlen/60, DT_hhmmss, DTend_hhmmss, 1.) 
 
  ax.set_title(title, size='x-large')
  ax.set_aspect('equal')

  at = AnchoredText("Max Probability: %d %s" % (prob1.max(),"%"), loc=4, prop=dict(size=10), frameon=True,)
  at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
  ax.add_artist(at)

  at = AnchoredText("%s" % (plotfilename), loc=2, prop=dict(size=8), frameon=True,)
  at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
  ax.add_artist(at)

# End Plot Section
#-----------------
# Write file

if _write_prob and not probfile:
  if threshold2 > 0.0 and xobs2.size > 0:
    write_prob(plotfilename, lat2d, lon2d, prob1, threshold, 1000., coards, starttime, endtime, \
               thres2 = threshold2, xlats = xobs2, ylats = yobs2)
  else:
    write_prob(plotfilename, lat2d, lon2d, prob1, threshold, 1000., coards, starttime, endtime)

# Save fig
#---------

if _yes_plot:
  P.savefig(plotfilename+".pdf",dpi=300)

  if not options.noshow:
    P.show()

#########################################################################
# Fortran90 Code for Cressman routine
#
# compile the fortran cressman.f90 routine with...
# 
# f2py --fcompiler="gnu95" --f90flags='-m64 -O3 -funroll-loops' -DF2PY_REPORT_ON_ARRAY_COPY -c -m cressman cressman.f90

#subroutine cressman(x, y, ob, xg, yg, roi, anal, nobs, nx, ny)
#  implicit none
#  integer :: nobs, nx, ny
#  real(8) :: x(nobs), y(nobs), ob(nobs)
#  real(8) :: xg(nx), yg(ny)
#  real(8), intent(out) :: anal(nx,ny)
#  real(8) :: roi
#
#!f2py real(8), intent(in),  dimension(nobs) :: x
#!f2py real(8), intent(in),  dimension(nobs) :: y
#!f2py real(8), intent(in),  dimension(nobs) :: ob
#!f2py real(8), intent(in),  dimension(nx)   :: xg
#!f2py real(8), intent(in),  dimension(ny)   :: yg
#!f2py real(8), intent(out), dimension(nx,ny) :: anal
#!f2py real,    intent(in) :: roi
#!f2py integer, intent(in) :: nobs, nx, ny
#
#  integer n, i, j
#  real dis, R2, w_sum, top, wk, rk2
#  real, parameter :: hsp = 1.33
#
#  R2 = roi**2.0
#
#! print *, 'Maxval of anal before:  ', maxval(anal)
#
#  DO j = 1,ny
#   DO i = 1,nx
#    w_sum = 0.0
#    top   = 0.0  
#!   anal(i,j) = 0.0
#    DO n = 1,nobs
#      dis = sqrt( (xg(i) - x(n))**2 + (yg(j)-y(n))**2 )
#      IF (dis .le. roi) THEN
#        rk2 = dis**2.0
#        wk = (R2-rk2) / (R2+rk2)
#        top = top + wk*ob(n)
#        w_sum = w_sum + wk
#      ENDIF
#
#!     IF (dis .le. 4.*roi) THEN
#!       wk = exp( -((dis/(hsp*roi))**2) )
#!       top = top + wk*ob(n)
#!       w_sum = w_sum + wk
#!     ENDIF
#
#    ENDDO
#
#    IF (w_sum .ge. 0.0001) THEN
#     anal(i,j) = anal(i,j) + min(top/w_sum,1.0)
#    ENDIF
#
#   ENDDO
#  ENDDO
#
#! print *, 'Maxval of anal after:  ', maxval(anal)
#end subroutine cressman
