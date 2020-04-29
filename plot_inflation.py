#!/usr/bin/env python
#
#
import matplotlib
import pylab as P
import numpy as N
import sys
import netCDF4
from optparse import OptionParser
from netCDF4 import num2date
import os
import ctables
import datetime as DT
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
import time as timeit
from cbook2 import *
import glob, pickle
from Plotting import shapefiles, viridis

interactive = True
output_format = "png"

# default time and z-level for plotting

# Other plotting stuff....
_inflation_ctable = viridis
_inflation_clevels = N.linspace(1.01, 5.01, num=21, endpoint=True)
_inflation_levels = [2, 10, 30]
_fig_size = (12,12)

#===============================================================================
def mymap(x, y, glat, glon, scale = 1.0, ax = None, ticks = True, resolution='c',\
          area_thresh = 10., shape_env = False, states=False, counties=False, pickle = False):

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

    if pickle:
        pickle.dump(map,open('mymap.pickle','wb'),-1)
        print(timeit.clock()-tt,' secs to create original Basemap instance and pickle it')

    return map
#===============================================================================
def plot_inflation(inflation, x, y, z, glat, glon, shape_env=None):

    fig, ((ax1, ax2), (ax3, ax4)) = P.subplots(2, 2, sharex=True, sharey=True, figsize = _fig_size)

    map = mymap(x, y, glat, glon, ax = ax1, shape_env=shape_env)

# get coordinates for contour plots

    lon2d, lat2d, xx, yy = map.makegrid(x.size, y.size, returnxy=True)

    nlevel = _inflation_levels[0]
    plot    = map.contourf(xx, yy, infl[nlevel], _inflation_clevels, cmap=_inflation_ctable)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    cbar.set_label("Inflation")
    plot    = map.contour(xx, yy,  infl[nlevel], _inflation_clevels[::2], colors='k', linewidths=0.5)
    title   = ("Inflation Factor Z = %3.2f km" % (z[nlevel]/1000.))
    ax1.set_title(title, fontsize=10)
#   if zoom:
#     ax1.set_xlim(1000*zoom[0],1000*zoom[1])
#     ax1.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max Inflation: %4.1f" % (infl[nlevel].max()), loc=4, prop=dict(size=6), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax1.add_artist(at)

    map = mymap(x, y, glat, glon, scale = hscale, ax = ax2, shape_env=shape_env)

# get coordinates for contour plots

    lon2d, lat2d, xx, yy = map.makegrid(x.size, y.size, returnxy=True)

    nlevel = _inflation_levels[1]
    plot    = map.contourf(xx, yy, infl[nlevel], _inflation_clevels, cmap=_inflation_ctable)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    cbar.set_label("Inflation")
    plot    = map.contour(xx, yy,  infl[nlevel], _inflation_clevels[::2], colors='k', linewidths=0.5)
    title   = ("Inflation Factor Z = %3.2f km" % (z[nlevel]/1000.))
    ax2.set_title(title, fontsize=10)
#   if zoom:
#     ax2.set_xlim(1000*zoom[0],1000*zoom[1])
#     ax2.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max Inflation: %4.1f" % (infl[nlevel].max()), loc=4, prop=dict(size=6), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax2.add_artist(at)

    map = mymap(x, y, glat, glon, scale = hscale, ax = ax3, shape_env=shape_env)

# get coordinates for contour plots

    lon2d, lat2d, xx, yy = map.makegrid(x.size, y.size, returnxy=True)

    nlevel = _inflation_levels[2]
    plot    = map.contourf(xx, yy, infl[nlevel], _inflation_clevels, cmap=_inflation_ctable)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    cbar.set_label("Inflation")
    plot    = map.contour(xx, yy,  infl[nlevel], _inflation_clevels[::2], colors='k', linewidths=0.5)
    title   = ("Inflation Factor Z = %3.2f km" % (z[nlevel]/1000.))
    ax3.set_title(title, fontsize=10)
#   if zoom:
#     ax2.set_xlim(1000*zoom[0],1000*zoom[1])
#     ax2.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max Inflation: %4.1f" % (infl[20].max()), loc=4, prop=dict(size=6), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax3.add_artist(at)

    map = mymap(x, y, glat, glon, scale = hscale, ax = ax4, shape_env=shape_env)

# get coordinates for contour plots

    lon2d, lat2d, xx, yy = map.makegrid(x.size, y.size, returnxy=True)

    plot    = map.contourf(xx, yy, infl.max(axis=0), _inflation_clevels, cmap=_inflation_ctable)
    cbar    = map.colorbar(plot,location='right',pad="5%")
    cbar.set_label("Inflation")
    plot    = map.contour(xx, yy,  infl.max(axis=0), _inflation_clevels[::2], colors='k', linewidths=0.5)
    title   = ("Max Inflation Factor in Column")
    ax4.set_title(title, fontsize=10)
#   if zoom:
#     ax2.set_xlim(1000*zoom[0],1000*zoom[1])
#     ax2.set_ylim(1000*zoom[2],1000*zoom[3])

    at = AnchoredText("Max Inflation: %4.1f" % (infl.max()), loc=4, prop=dict(size=6), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax4.add_artist(at)
    
    return fig

#---------------------------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls
#
parser = OptionParser()

parser.add_option("-e",  "--exp", dest="exp", type="string", help = "Path to the pickled experiment file generated \
                                                                     from the create_run_letkf script")
parser.add_option("-d",  "--dir", dest="dir", type="string", help = "Experiment directory - assumes that the experiment \
                                                                     database file is named the same as directory")
parser.add_option("-t",  "--time",   dest="datetime", default=None, type = "string", help = "Usage:  --time 2003,5,8,21,0,0")
parser.add_option("--show", dest="show", action="store_true", help="Show plot interactively")

(options, args) = parser.parse_args()

if options.dir and options.exp == None:
    options.exp = glob.glob(os.path.join(options.dir, "*.exp" ))[0]
    print("\n ==> Plot_Inflation: found experiment files %s" % options.exp)

if options.exp == None:
    parser.print_help()
    print("\n ==> Plot_Inflation: ERROR --> Experiment file not supplied..EXITING!!!")
    sys.exit(-1)
else:
    with open(options.exp, 'rb') as f:
        exper = pickle.load(f)

if options.datetime == None:
    parser.print_help()
    print("\n ==> Plot_Inflation: ERROR --> time for files not supplied..EXITING!!!")
    sys.exit(-1)
else:
    list = [int(t) for t in options.datetime.split(",")]
    myDT = datetime.datetime(list[0],list[1],list[2],list[3],list[4],list[5])
    print("\n ==> Plot_Inflation: Date and time supplied is %s" % (myDT.strftime("%Y-%m-%d %H:%M:%S")))

infilename = "Inflation_%s.nc" % (myDT.strftime("%Y-%m-%d_%H:%M:%S"))
newfile = os.path.join(exper['base_path'], infilename)  

print("\n ==> Plot_Inflation: opening inflation file:  %s" % (newfile))

file_obj  = netCDF4.Dataset(newfile, "r")
xc        = file_obj.variables['XC'][:]
yc        = file_obj.variables['YC'][:]
zc        = file_obj.variables['ZC'][:]
infl      = file_obj.variables['inflation'][...]
glat = exper['lat0']
glon = exper['lon0']

fig = plot_inflation(infl, xc, yc, zc, glat, glon, shape_env=shapefiles)

title = ("Inflation at %s " % (myDT.strftime("%Y-%m-%d  %H:%M:%S")))
fig.suptitle(title, fontsize=16)

# Write out file into the "Plots" directory
newbase = os.path.join("./", os.path.split(newfile)[0], "Plots", os.path.split(newfile)[1][:-2])

filename = "%s%s" % (newbase, "pdf")
print "\n Saving file %s" % (filename)
fig.savefig(filename, format="pdf", dpi=300)
 
if options.show:
   os.system("open %s" % filename)
 
# End of file
