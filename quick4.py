#!/usr/bin/env python
#
#
import matplotlib
import pylab as P
import numpy as N
import sys
import glob
import pickle
import netCDF4 as ncdf
from optparse import OptionParser
import os
import ctables
import datetime as DT
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
import time as timeit
from cbook2 import *
from Plotting import shapefiles

_debug_timing = False

# Plotting defaults
output_format = "pdf"
_fig_size = (12,12)

# z-level for plotting
_height  = 1000.
#default member
_member  = 1

_min_dbz = 10.
_min_w   = 0.1

# Some color tables and scales for vorticity...
_ref_ctable = ctables.REF_default
_wz_clevels = N.arange(-150.,175.,25.)
_w_clevels  = N.arange(-15.,16.,1.)

#=======================================================================================================================
# Find the CM1 file for time requested
def FindRestartFile(fcst_path,run_name,time):

    fileheader = os.path.join(fcst_path,run_name)
    files = glob.glob(fileheader+"_rst_0*.nc")

    if len(files) > 0:

        for n, file in enumerate(files):
            f = ncdf.Dataset(file, "r")
            f_time = f.variables['time'][0]

            if( N.abs(f_time - time) < 1.0 ):
                print("\n Found time %d in file:  %s" % (f_time, file))
                return n

        print("\n\n  ERROR!!!")
        print("\n      FindRestartFile: A file with the restart time of %d cannot be found, exiting!!!" % time)
        print("\n\n  ERROR!!!")
        sys.exit(-1)

    else:
        print("\n\n  ERROR!!!")
        print("\n      FindRestartFile: No files where found having the header:  %s, exiting" % fileheader)
        print("\n\n  ERROR!!!")
        sys.exit(-1)


#===============================================================================
def mymap(x, y, glat, glon, scale = 1.0, ax = None, ticks = True, resolution='c',\
          area_thresh = 10., shape_env = False, states=False, counties=False, pickle = False):

    if _debug_timing:
        tt = timeit.clock()

    xmax = max(x) / scale
    xmin = min(x) / scale
    ymax = max(y) / scale
    ymin = min(y) / scale

    sw_lat, sw_lon = dxy_2_dll(xmin, ymin, glat, glon, degrees=True)
    ne_lat, ne_lon = dxy_2_dll(xmax, ymax, glat, glon, degrees=True)
        
    map = Basemap(llcrnrlon=sw_lon, llcrnrlat=sw_lat, \
                  urcrnrlon=ne_lon, urcrnrlat=ne_lat, \
                  lat_0=0.5*(ne_lat+sw_lat), lon_0=0.5*(ne_lon+sw_lon), \
                  projection = 'lcc',      \
                  resolution=resolution,   \
                  area_thresh=area_thresh, \
                  suppress_ticks=ticks, \
                  ax=ax)

    if counties:
        map.drawcounties()
        
    if states:
        map.drawstates()


# Shape file stuff

    if shape_env:

        for item in shape_env:
            items = item.split(",")
            shapefile  = items[0]
            color      = items[1]
            linewidth  = items[2]

            s = map.readshapefile(shapefile,'shapeinfo',drawbounds=False)

            for shape in map.shapeinfo:
                xx, yy = zip(*shape)
                map.plot(xx,yy,color=color,linewidth=linewidth,ax=ax)

# pickle the class instance.

    if _debug_timing:
        print(timeit.clock()-tt,' secs to create original Basemap instance')

    if pickle:
        pickle.dump(map,open('mymap.pickle','wb'),-1)
        print(timeit.clock()-tt,' secs to create original Basemap instance and pickle it')

    return map
#===============================================================================
def plot_W_DBZ_T_WZ(w, dbz, t, wz, x, y, height, time, member, glat=None, glon=None, sfc=False, \
                    filename=None, nodisplay=False, shape_env = False, zoom=None):

    if filename == None:
        filename = "%s_%4.2f_MEMBER_%2.2d.%s" % (time.strftime("%Y_%m_%d_%H:%M:%S"), height, member, output_format)
    else:
        filename = "%s_%s_%4.2f_MEMBER_%2.2d.%s" % (filename, time.strftime("%H:%M:%S"), height, member)
               
    fig, ((ax1, ax2), (ax3, ax4)) = P.subplots(2, 2, sharex=True, sharey=True, figsize = _fig_size)

    map = mymap(x, y, glat, glon, ax = ax1, shape_env=shape_env)

# get coordinates for contour plots

    lon2d, lat2d, xx, yy = map.makegrid(x.size, y.size, returnxy=True)

    clevels = N.arange(0.,75.,5.)
    plot    = map.contourf(xx, yy, N.ma.masked_less_equal(dbz,_min_dbz), clevels, cmap=_ref_ctable)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    cbar.set_label("dBZ")
    plot    = map.contour(xx, yy,  dbz, clevels[::2], colors='k', linewidths=0.5)
    title   = ("Reflectivity")
    ax1.set_title(title, fontsize=12)
    if zoom:
      ax1.set_xlim(1000*zoom[0],1000*zoom[1])
      ax1.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max dBZ: %4.1f" % (dbz.max()), loc=4, prop=dict(size=12), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax1.add_artist(at)

    map = mymap(x, y, glat, glon, scale = hscale, ax = ax2, shape_env=shape_env)

    scale_w_clevels = min(max(N.int(height/1000.), 1.0), 7.0)
    clevels = scale_w_clevels*N.arange(-15.,16.,1.)
    wmask   = N.ma.masked_array(w, mask = [N.abs(w) <= scale_w_clevels*_min_w])
    plot    = map.contourf(xx, yy, wmask, clevels, cmap=ctables.Not_PosDef_Default)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    plot    = map.contour(xx, yy, wmask, clevels[::2], colors='k', linewidths=0.5)
    cbar.set_label('%s' % ("$m s^{-1}$"))
    title = ("Vertical Velocity")
    ax2.set_title(title, fontsize=12)
    if zoom:
      ax2.set_xlim(1000*zoom[0],1000*zoom[1])
      ax2.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max W: %4.1f \n Min W: %4.1f" % (w.max(),w.min()), loc=4, prop=dict(size=12), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax2.add_artist(at)

    map = mymap(x, y, glat, glon, scale = hscale, ax = ax3, shape_env=shape_env)

    clevels = N.arange(-10.,11.,1.)
    plot    = map.contourf(xx, yy, t, clevels, cmap=ctables.Not_PosDef_Default)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    plot    = map.contour(xx, yy, t, clevels[::2], colors='k', linewidths=0.5)
    cbar.set_label('%s' % ("K"))
    if sfc:
        title = ("SFC Pert. Potential Temperature")
    else:
        title = ("Pert. Potential Temperature")
    ax3.set_title(title, fontsize=12)
    if zoom:
      ax3.set_xlim(1000*zoom[0],1000*zoom[1])
      ax3.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max TH: %4.1f \n Min TH: %4.1f" % (t.max(),t.min()), loc=4, prop=dict(size=12), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax3.add_artist(at)

    map = mymap(x, y, glat, glon, scale = hscale, ax = ax4, shape_env=shape_env)

    s_wz    = wz*10000.
    plot    = map.contourf(xx, yy, s_wz, _wz_clevels, cmap=ctables.Not_PosDef_Default)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    plot    = map.contour(xx, yy, s_wz, _wz_clevels[::2], colors='k', linewidths=0.5)
    cbar.set_label('%s' % ("x $ 10^{4}s^{-1}$"))
    if sfc:
        title = ("SFC Vert. Vorticity")
    else:
        title = ("Vert. Vorticity")
    ax4.set_title(title, fontsize=12)
    if zoom:
      ax4.set_xlim(1000*zoom[0],1000*zoom[1])
      ax4.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max Wz: %4.1f \n Min Wz: %4.1f" % (s_wz.max(),s_wz.min()), loc=4, prop=dict(size=12), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax4.add_artist(at)

    title = ("\n %s  Height:  %4.2f km      MEMBER = %d" % \
            (time.strftime("%Y_%m_%d  %H:%M:%S"),height/1000.,member))
    fig.suptitle(title, fontsize=16)

    return fig, filename

#---------------------------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls
#

parser = OptionParser()

parser.add_option("-e",  "--exp", dest="exp", type="string", help = "Path to the pickled experiment file generated \
                                                                     from the create_run_letkf script")
parser.add_option("-d",  "--dir", dest="dir", type="string", help = "Experiment directory - assumes that the experiment \
                                                                     database file is named the same as directory")

parser.add_option("-t",  "--time",   dest="datetime", default=None, type = "string", help = "Usage:  --time 2003,5,8,21,0,0")

parser.add_option("-o", "--output", dest="outfile", type="string", default= None, \
                                    help="Name of file to be written out (no suffix included)")
parser.add_option("-z", "--height", dest="height", type="float", help = "Model height in meters...1000., 2000., ")

parser.add_option("-s", "--sfc",    dest="sfc",    action="store_true", help="plot surface vorticity and pert pot temp")

parser.add_option("--nodisplay", dest="nodisplay", action="store_true", help="dont show plot interactively")

parser.add_option(       "--zoom",     dest="zoom", type="int", nargs = 4, help="bounds (km) of plot - 4 args required: xmin xmax ymin ymax)")

parser.add_option("-m", "--member",    dest="member", type="int" )

(options, args) = parser.parse_args()

if options.dir and options.exp == None:
    options.exp = glob.glob(os.path.join(options.dir, "*.exp" ))[0]
    print("\n ==> Quick4Panel: found experiment files %s" % options.exp)

if options.exp == None:
    parser.print_help()
    print("\n ==> Quick4Panel: ERROR --> Experiment file not supplied..EXITING!!!")
    sys.exit(-1)
else:
    with open(options.exp, 'rb') as f:
        exper = pickle.load(f)

if options.datetime == None:
    parser.print_help()
    print("\n ==> Quick4Panel: ERROR --> time for files not supplied..EXITING!!!")
    sys.exit(-1)
else:
    list = [int(t) for t in options.datetime.split(",")]
    myDT = datetime.datetime(list[0],list[1],list[2],list[3],list[4],list[5])
    print("\n ==> Quick4Panel: Date and time supplied is %s" % (myDT.strftime("%Y %m-%d %H:%M:%S")))

if options.height == None:
    print
    print("\n ==> Quick4Panel:  No height supplied, using default height: %4.2f\n" % _height)
    print
    height = _height
else:
    height = options.height

if options.member == None:
    print
    print("\n ==> Quick4Panel:  No member supplied, using default member: %d\n" % _member)
    print
    member = _member
else:
    member = options.member-1

# Figure out file time

mDT  = datetime.datetime(exper['YEAR'],
                         exper['MONTH'],
                         exper['DAY'],
                         exper['HOUR'],
                         exper['MINUTE'],
                         exper['SECOND'])

time = (myDT - mDT).seconds
number   = FindRestartFile(exper['fcst_members'][0],"cm1out",time)

files = []

for f in exper['fcst_members']:
    newBasePath = os.path.join("./", os.path.split(options.exp)[0])
    files.append(os.path.join(newBasePath, os.path.split(f)[1], ("cm1out_rst_%6.6d.nc" % number)))

options.file = files[member]

f         = netCDF4.Dataset(options.file, "r")
glat      = exper['lat0']
glon      = exper['lon0']
xoffset   = exper['xoffset']
yoffset   = exper['yoffset']

xc        = f.variables['xh'][:] + xoffset
yc        = f.variables['yh'][:] + yoffset
zc        = f.variables['zh'][:]
ze        = f.variables['zf'][:]

zb, zt, dzb, dzt, dz = interp_weights(height, ze)

wplot     = (f.variables['wa'][zb]*dzb + f.variables['wa'][zt]*dzt) / dz
zb, zt, dzb, dzt, dz = interp_weights(height, zc)
dplot     = (f.variables['dbz'][zb]*dzb + f.variables['dbz'][zt]*dzt) / dz

if options.sfc:
    tplot     = f.variables['theta'][0] - f.variables['th0'][0]
    wzplot    = ComputeWZ(xc, yc, f.variables['ua'][0], f.variables['va'][0])
else:
    tplot     = ((f.variables['theta'][zb]*dzb + f.variables['theta'][zt]*dzt) / dz) \
              - ((f.variables['th0'][zb]*dzb + f.variables['th0'][zt]*dzt) / dz)
    wzplot    = (ComputeWZ(xc, yc, f.variables['ua'][zb], f.variables['va'][zb])*dzb 
                +ComputeWZ(xc, yc, f.variables['ua'][zt], f.variables['va'][zt])*dzt) / dz

figure, filename = plot_W_DBZ_T_WZ(wplot, dplot, tplot, wzplot, xc, yc, height, myDT, \
                                   member=member, glat=glat, glon=glon, sfc=options.sfc, filename=options.outfile, \
                                   nodisplay=options.nodisplay, zoom=options.zoom,shape_env=shapefiles)

# Store plots in RUN EXP directory

newfilename = os.path.join(newBasePath,'Plots',filename)
if output_format: 
    print "\n Saving file %s" % (newfilename)
    figure.savefig(newfilename, format=output_format)

if not options.nodisplay:
    os.system("open %s" % newfilename)

# End of file
