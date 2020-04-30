#!/usr/bin/env python

import numpy as N
import netCDF4 as ncdf
import datetime as DT 
from optparse import OptionParser
import sys, os
import random
from assim_util import *
from scipy.stats import binned_statistic_2d
from scipy.interpolate import griddata
import pylab as P

# 
_missing_value = 9900000000.
_time_units    = 'seconds since 1970-01-01 00:00:00'
_calendar      = 'standard'
_stringlen     = 8
_datelen       = 19

#####################################################################################################
def create_prior_file(ne):
        
# create the fileput filename and create new netCDF4 file

  filename = "%s%s" % ("Prior", ".nc" )
    
  rootgroup = ncdf.Dataset(filename, 'w', format='NETCDF4')
      
# Create dimensions

  rootgroup.createDimension('ob_num', None)
  rootgroup.createDimension('stringlen', _stringlen)
  rootgroup.createDimension('datelen', _datelen)
  rootgroup.createDimension('NE', ne)
  
# Write some attributes

  rootgroup.time_units = _time_units
  rootgroup.calendar   = _calendar
  rootgroup.stringlen  = "%d" % (_stringlen)
  rootgroup.datelen    = "%d" % (_datelen)

# Create variables

  V_type   = rootgroup.createVariable('type', 'S1', ('ob_num', 'stringlen'), zlib=True, shuffle=True )
  
# V_units  = rootgroup.createVariable('units', 'S1', ('ob_num'), zlib=True, shuffle=True)
  
  V_secs   = rootgroup.createVariable('secs', 'i8', ('ob_num',), zlib=True, shuffle=True)
  V_secs.units = _time_units
  
  V_dates  = rootgroup.createVariable('dates', 'S1', ('ob_num', 'datelen'), zlib=True, shuffle=True)

#  
  V_kind     = rootgroup.createVariable('kind',  'i4', ('ob_num',))

#
  V_year     = rootgroup.createVariable('year',  'i4', ('ob_num',))
  V_month    = rootgroup.createVariable('month', 'i4', ('ob_num',))
  V_day      = rootgroup.createVariable('day',   'i4', ('ob_num',))
  V_hour     = rootgroup.createVariable('hour',  'i4', ('ob_num',))
  V_minute   = rootgroup.createVariable('minute', 'i4', ('ob_num',))
  V_second   = rootgroup.createVariable('second', 'i4', ('ob_num',))
      
  V_value = rootgroup.createVariable('value', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_value._missing_value = _missing_value 

  V_value = rootgroup.createVariable('error', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  
  V_x = rootgroup.createVariable('x', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_x._missing_value = _missing_value 
  V_x.units = "m"
  
  V_y = rootgroup.createVariable('y', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_y._missing_value = _missing_value 
  V_y.units = "m"
  
  V_z = rootgroup.createVariable('z', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_z._missing_value = _missing_value 
  V_z.units = "m"
  
  V_lat = rootgroup.createVariable('lat', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_lat._missing_value = _missing_value 
  V_lat.units = "deg"
  
  V_lon = rootgroup.createVariable('lon', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_lon._missing_value = _missing_value 
  V_lon.units = "deg"
  
  V_az = rootgroup.createVariable('az', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_az._missing_value = _missing_value 
  V_az.units = "deg"
  
  V_el = rootgroup.createVariable('el', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_el._missing_value = _missing_value 
  V_el.units = "deg"
  
  V_Hxf = rootgroup.createVariable('Hxf', 'f8', ('ob_num','NE'), zlib=True, shuffle=True) 
  V_Hxf._missing_value = _missing_value 
  
  V_Hxfbar = rootgroup.createVariable('Hxfbar', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_Hxfbar._missing_value = _missing_value 
  
  V_status   = rootgroup.createVariable('status', 'i4', ('ob_num',), zlib=True, shuffle=True)
  V_status._missing_value = _missing_value 
  
  V_value = rootgroup.createVariable('outlier', 'f8', ('ob_num',), zlib=True, shuffle=True) 
  V_value._missing_value = _missing_value 

  rootgroup.close()
##################################################################################################### 
def write_prior(ne, kind, value, dates, error, x, y, z, Hxf, Hxfbar, lat, lon, el, az, status, outlier):
  
  filename = "%s%s" % ("Prior", ".nc" )
    
  if os.path.isfile(filename) != True:
      print("No Prior file found, creating %s \n" % filename)
      create_prior_file(ne)
  else:
      print("Prior obs file found, opening %s \n" % filename)
  
# create the fileput filename and create new netCDF4 file

  rootgroup = ncdf.Dataset(filename, 'a', format='NETCDF4')
  
# get current number of obs

  oo = len(rootgroup.dimensions['ob_num'])
  stringlen = len(rootgroup.dimensions['stringlen'])
  datelen   = len(rootgroup.dimensions['datelen'])
  
# get number of obs to add 

  obs = N.size(value)
  
  print("Previous number of obs %d \n" % oo)
  print("Adding %d to netcdf file \n" % obs)
  
# get kind then setup type and units

  ltypes = N.empty([obs,stringlen], dtype=N.str)

  for k, ktype in enumerate(kind):
      if ktype == 12: 
          ltypes[k,:] = ncdf.stringtoarr("REFL", stringlen)
      if ktype == 11: 
          ltypes[k,:] = ncdf.stringtoarr("VR", stringlen)
      if ktype == 111: 
          ltypes[k,:] = ncdf.stringtoarr("dVR", stringlen)
          
  rootgroup.variables['type'][oo:,:] = ltypes
  rootgroup.variables['kind'][oo:]   = kind
      
# Write a date_time string AND separate the date_time info out individually.
# the table stores date as YYYY-MM-DD HH:mm:SS for example "2007-01-20 01:35:00"

  str_dates = []
  for n, itemb in enumerate(dates):
    item = itemb.decode("utf-8")
    rootgroup.variables['dates'][oo+n,:] = ncdf.stringtoarr(item, datelen)
    str_dates.append(str(item))
    try:
      dt = DT.datetime.strptime(str(item), "%Y-%m-%d_%H:%M:%S")
    except:
      dt = DT.datetime.strptime(str(item), "%Y-%m-%d %H:%M:%S")
    rootgroup.variables['secs'][oo+n] = N.array(ncdf.date2num(dt,units=_time_units,calendar=_calendar))

  rootgroup.variables['year'][oo:]   = [date[0:4]   for date in str_dates]
  rootgroup.variables['month'][oo:]  = [date[5:7]   for date in str_dates]
  rootgroup.variables['day'][oo:]    = [date[8:10]  for date in str_dates]
  rootgroup.variables['hour'][oo:]   = [date[11:13] for date in str_dates]
  rootgroup.variables['minute'][oo:] = [date[14:16] for date in str_dates]
  rootgroup.variables['second'][oo:] = [date[17:19] for date in str_dates]

  rootgroup.variables['value'][oo:]  = value
  rootgroup.variables['error'][oo:]  = error   # this is the variance, not std. deviation!!!
  
  rootgroup.variables['x'][oo:] = x
  rootgroup.variables['y'][oo:] = y
  rootgroup.variables['z'][oo:] = z
  rootgroup.variables['lat'][oo:] = lat
  rootgroup.variables['lon'][oo:] = lon
  rootgroup.variables['el'][oo:] = el
  rootgroup.variables['az'][oo:] = az
  rootgroup.variables['status'][oo:] = status
  rootgroup.variables['outlier'][oo:] = outlier
          
  rootgroup.variables['Hxf'][oo:,:]  = Hxf[:,:]
  rootgroup.variables['Hxfbar'][oo:] = Hxfbar
 
  rootgroup.sync()
  
  print("Total number of obs in %s file is %d \n" % (filename,len(rootgroup.dimensions['ob_num'])))
  rootgroup.close()

##################################################################################################### 
def read_prior(file_DT, file=None, sort=False, dict=False):
  
  if file != None:
    filename = file
  else:
    filename = "%s_%s%s" % ("Prior", file_DT.strftime("%Y-%m-%d_%H:%M:%S"), ".nc")
    
  if os.path.isfile(filename) != True:
      print("\n==> No Prior file found:  %s ... exiting!!! \n" % filename)
      sys.exit(-1)
  else:
      print("Prior obs file found, opening %s \n" % filename)
  
# create the fileput filename and create new netCDF4 file

  rootgroup = ncdf.Dataset(filename, 'r')
  
# get current number of obs

  oo = len(rootgroup.dimensions['ob_num'])
  
  value = rootgroup.variables['value'][:]  
  error = rootgroup.variables['error'][:]  # this is the variance, not std. deviation!!!
  
  x     = rootgroup.variables['x'][:] 
  y     = rootgroup.variables['y'][:] 
  z     = rootgroup.variables['z'][:] 
  lat   = rootgroup.variables['lat'][:] 
  lon   = rootgroup.variables['lon'][:]
  el    = rootgroup.variables['el'][:] 
  az    = rootgroup.variables['az'][:] 
  status= rootgroup.variables['status'][:]
  outlier= rootgroup.variables['outlier'][:]
          
  Hxf    = rootgroup.variables['Hxf'][:,:]  
  Hxfbar = rootgroup.variables['Hxfbar'][:] 
  secs   = rootgroup.variables['secs'][:]
  kind   = rootgroup.variables['kind'][:]
  dates  = rootgroup.variables['dates'][...]
 
  print("Total number of obs in %s file is %d \n" % (filename,oo))
  
  rootgroup.close()
  
# Create a temporay record array to use the numpy sort capability - then use the indices generated in 'index'
#        to return a the standard set of data sorted by x/y/z/sec/kind to the calling routine

  if sort:
    array = N.recarray((x.size,), dtype = [('x', N.float64), ('y', N.float64), ('z', N.float64), \
                                           ('lat', N.float64), ('lon', N.float64), ('secs', N.int64), \
                                           ('kind', N.int),  ('index', N.int64), ('outlier',N.float64)])
    array['x']       = x
    array['y']       = y
    array['z']       = z
    array['lat']     = lat
    array['lon']     = lon
    array['secs']    = secs
    array['kind']    = kind
    array['outlier'] = outlier
    array['index']   = N.arange(x.size)

    array.sort(order=['x','y','z','secs','kind'])

    return kind[array['index']], value[array['index']], error[array['index']], \
           lat[array['index']], lon[array['index']], x[array['index']], y[array['index']], z[array['index']], \
           secs[array['index']], Hxf[array['index'],:], Hxfbar[array['index']], outlier[array['index']]

  else:
    if dict:
      return {'kind': kind, 'value': value, 'error': error,
              'x': x,        'y': y,         'z': z,     'el': el,          'az': az,
              'lat': lat,    'lon': lon,    'secs': secs,  'status':  status,
              'Hxf': Hxf, 'Hxfbar': Hxfbar, 'dates': ncdf.chartostring(dates), 'outlier': outlier}
    else:
      return kind, value, error, lat, lon, x, y, z, secs, Hxf, Hxfbar, outlier

#####################################################################################################
def read_prior_for_testing(analysis_time, prior_file):

  outlier_threshold = 3

  analysis_time = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))

  kind, value, error, lat, lon, xob, yob, zob, tob, Hxf, Hxfbar, outlier = read_prior(analysis_time, file=prior_file, sort=True)

  nobs = Hxf.shape[0]
  ne   = Hxf.shape[1]

  if nobs == 0:
    print("---->Prior file is empty, exiting at analysis time %s \n" %  (analysis_time.strftime("%Y %m-%d %H:%M:%S")))
    sys.exit(0)
  else:
    print("---->Prior file has %d obs at analysis time %s \n" %  (nobs, analysis_time.strftime("%Y %m-%d %H:%M:%S")))

    print('n\test called with an outlier threshold of %d standard deviations' % outlier_threshold)

    mask = N.where( outlier <= outlier_threshold, True, False )
    print(nobs, N.count_nonzero(mask))
    mask2 = N.where( kind == 11, True, False )
    mask  = mask | mask2
    print(nobs, N.count_nonzero(mask2), N.count_nonzero(mask))

#####################################################################################################
def compute_dbz_bias(time, prior_file, dbz_threshold=10., plot=False):

  outlier_threshold = 3

  analysis_time = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))

  kind, value, error, lat, lon, xob, yob, zob, tob, Hxf, Hxfbar, outlier = read_prior(analysis_time, file=prior_file, sort=True)

  nobs = Hxf.shape[0]
  ne   = Hxf.shape[1]

  if nobs == 0:
    print("---->Prior file is empty, exiting at analysis time %s \n" %  (analysis_time.strftime("%Y %m-%d %H:%M:%S")))
    sys.exit(0)
  else:
    print("---->Prior file has %d obs at analysis time %s \n" %  (nobs, analysis_time.strftime("%Y %m-%d %H:%M:%S")))


#  BIAS CORRECTION
#  FIND OUT WHERE DBZ exists and is greater than some threshold
    mask = N.where( kind == 12, True, False )
    print("  --> COMPUTE_DBZ_BIAS:  Variable dBZ error:  Number of dBZ observations is %d" % N.count_nonzero(mask))
    mask2 = N.where(value > dbz_threshold, True, False)
    print("  --> COMPUTE_DBZ_BIAS:  Variable dBZ error:  Number of obs where value > %d  %d" % (dbz_threshold, N.count_nonzero(mask2)))
    mask  = mask & mask2
    print("  -->  COMPUTE_DBZ_BIAS:  Number of Priors where Obs > %d and will be used for bias correction:  %d" % (dbz_threshold, N.count_nonzero(mask)))

# CREATE BINS FOR GRIDDING
    z_bins   = N.arange(0., 11000., 1000.)    # Z-bins
    dbz_bins = N.arange(dbz_threshold, 80., 10.)
    values   = value[mask] - Hxfbar[mask]
    Hxftmp   = Hxfbar[mask]
    
    if plot:
      fig = P.figure(figsize=(8,12))
      P.subplot(3, 1, 1)
      grid_x, grid_y = N.meshgrid(dbz_bins,z_bins)
    
      stats = griddata((Hxftmp,zob[mask]), values, (grid_x, grid_y), method='nearest')
      extent = [dbz_bins[0], dbz_bins[-1], 0.001*z_bins[0], 0.001*z_bins[-1]]
      im = P.imshow(stats.transpose(), cmap='RdBu_r', extent=extent, interpolation='nearest', origin='lower')
      P.colorbar(im, extend='both', orientation='horizontal', shrink=0.8)
      P.title("Raw interpolation using grid data")

      P.subplot(3, 1, 2)
      stats, xedges, yedges, bins = binned_statistic_2d(Hxftmp, zob[mask], values, statistic='mean', bins=[dbz_bins,z_bins])
    
      extent = [xedges[0], xedges[-1], 0.001*yedges[0], 0.001*yedges[-1]]
      im = P.imshow(stats, cmap='RdBu_r', extent=extent, interpolation='nearest', origin='lower')
      P.colorbar(im, extend='both', orientation='horizontal', shrink=0.8)
      P.axis(extent)
      P.title("Bins using scipy.binned_statistic_2d")

      P.subplot(3, 1, 3)
      ne       = Hxf.shape[1]
      values   = N.repeat(value[mask], ne) - Hxf[mask,:].flatten()
      Hxftmp   = Hxf[mask,:].flatten()
      zobs     = N.repeat(zob[mask], ne)
      stats, xedges, yedges, bins = binned_statistic_2d(Hxftmp, zobs, values, statistic='mean', bins=[dbz_bins,z_bins])
    
      extent = [xedges[0], xedges[-1], 0.001*yedges[0], 0.001*yedges[-1]]
      im = P.imshow(stats, cmap='RdBu_r', extent=extent, interpolation='nearest', origin='lower')
      P.colorbar(im, extend='both', orientation='horizontal', shrink=0.8)
      P.axis(extent)
      P.title("Bins using scipy.binned_statistic_2d for ALL priors")

#   P.scatter(Hxftmp, 0.001*zob[mask], s=5, c='0.5', edgecolor='0.5')
      P.show()
    
    return binned_statistic_2d(Hxftmp, zobs, values, statistic='mean', bins=[dbz_bins,z_bins])
#####################################################################################################
# Main program
if __name__ == "__main__":

  now = DT.datetime.now()

  print("----------------------------------------------------------------------\n")
  print("              BEGIN PROGRAM PRIOR_FILE                             \n ")
  print("  WALLCLOCK START TIME:  %s \n" % now.strftime("%Y-%m-%d %H:%M")  )
  print("  --------------------------------------------------------------------\n")
  

# Parse input command options

  parser = OptionParser()
  parser.add_option(      "--dt",         dest="time", type="string",  help = "date and time")
  parser.add_option(      "--test",       dest="test", default=False,  help = "create fake data and rd/wt", action="store_true")
  parser.add_option(      "--test2",      dest="test2", default=False,  help = "call dummy test routine", action="store_true")
  parser.add_option(      "--bias",       dest="bias", default=False,  help = "Test Prior bias analysis", action="store_true")

  parser.add_option("-f", "--file",   dest="file",     type="string", help = "Priors netCDF file")
  parser.add_option("-t", "--time",   dest="time",     type="string", help = "Analysis time in YYYY,MM,DD,HH,MM,SS")

  (options, args) = parser.parse_args()
 
# TEST CODE FOR BIAS CORRECTION
 
  if options.bias:
    if options.time == None:
      print("\n ====>No analysis time specified, exiting....\n")
      parser.print_help()
      sys.exit(-1)
    else:
      time = options.time.split(",")

    if options.file == None:
      print("\n ====>No netCDF Priors file specified, using default priors file\n")
      ftime      = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))
      prior_file = "Prior_%s.nc" % ftime.strftime("%Y-%m-%d_%H:%M:%S")
      if os.path.isfile(prior_file):
        print("\n Using prior file:  %s\n" % prior_file)
      else:
        print("\n ====>Could not find prior file:  %s !!!!  Exiting LETKF....\n" % prior_file)
        parser.print_help()
        sys.exit(-1)
    else:
      prior_file = options.file

    ret = compute_dbz_bias(time, prior_file, plot=True)

#
# More test code

  if options.test2 == True:
    if options.time == None:
      print("\n ====>No analysis time specified, exiting....\n")
      parser.print_help()
      sys.exit(-1)
    else:
      time = options.time.split(",")

    if options.file == None:
      print("\n ====>No netCDF Priors file specified, using default priors file\n")
      ftime      = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))
      prior_file = "Prior_%s.nc" % ftime.strftime("%Y-%m-%d_%H:%M:%S")
      if os.path.isfile(prior_file):
        print("\n Using prior file:  %s\n" % prior_file)
      else:
        print("\n ====>Could not find prior file:  %s !!!!  Exiting LETKF....\n" % prior_file)
        parser.print_help()
        sys.exit(-1)
    else:
      prior_file = options.file
      

    read_prior_for_testing(time, prior_file)

  if options.test == True:
    print("\n Testing nc_prior_file\n")
    ne = 20
    options.time = now
  
# Create fake data....

    create_prior_file(ne)
  
    d = DT.datetime(2003,5,8,22,10,0)

    kind   = N.array((11,12), dtype = N.int32)
    status = N.array((1,1), dtype = N.int32)
    value  = N.array((-23.2,34.6), dtype = N.float32)
    dates  = N.array((d.strftime("%Y-%m-%d %H:%M:%S"),d.strftime("%Y-%m-%d %H:%M:%S")))
    error  = N.array((2.0,5.0), dtype = N.float32)
    x      = N.array((10000.,1000.), dtype = N.float)
    y      = N.array((-10000.,1000.), dtype = N.float)
    lat    = N.array((98.0,98.0), dtype = N.float)
    lon    = N.array((101.1,101.1), dtype = N.float)
    z      = N.array((320.,320.), dtype = N.float)
    az     = N.array((33.,33.), dtype = N.float)
    el     = N.array((0.5,0.5), dtype = N.float)
    Hxf    = N.zeros((2,ne), dtype = N.float32)
    Hxfbar = N.zeros((2), dtype = N.float32)
    outlier= N.array((4.0,0.5), dtype = N.float)
  
    for n in N.arange(ne):
      Hxf[0,n] = value[0] + error[0] * random.uniform(-0.5,0.5)
      Hxf[1,n] = value[1] + error[1] * random.uniform(-0.5,0.5)
    
    Hxfbar[0] = Hxf[0,:].mean()
    Hxfbar[1] = Hxf[1,:].mean()

# Write fake data out...

    write_prior(ne, kind, value, dates, error, x, y, z, Hxf, Hxfbar, lat, lon, el, az, status, outlier)
 
    newPriorFile = "%s_%s%s" % ("Prior", d.strftime("%Y-%m-%d_%H:%M:%S"), ".nc")
    os.rename("Prior.nc", newPriorFile)
    print("\n==> Moved Prior.nc file to %s\n" % (newPriorFile))
  
# read fake data in....

    kind0, value0, error0, lat0, lon0, x0, y0, z0, t0, Hxf0, Hxfbar0, outlier = read_prior(d)
  
    print("Input Hxf - Output Hxf == Zero?")
  
    print(Hxf - Hxf0)

    print("Kinds:  ", kind0)

