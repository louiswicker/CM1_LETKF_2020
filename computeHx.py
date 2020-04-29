#!/usr/bin/env python
import sys
import os
import numpy as N
import datetime as DT 
from netCDF4 import *
import pyDart
import ens 
import glob
from optparse import OptionParser
from multiprocessing import Pool
from time import time as timer
from nc_prior_file import *
from assim_util import *
import state_vector as state
import json

missing = -999.
diag = True

# Stuff to set observation std deviation in filter

_default_obs_error = {11: ["VR", 3.0], 12:["DBZ", 7.5]}
  
#-----------------------------------------------------------------------------------------------------------------------------------
# Main program

if __name__ == "__main__":

  stime = timer()
  now = DT.datetime.now()

  print("\n  ----------------------------------------------------------------------")
  print("\n                BEGIN PROGRAM ComputeHx                                 ")
  print("\n    WALLCLOCK START TIME:  %s \n" % now.strftime("%Y-%m-%d %H:%M")  )
  print("\n  ----------------------------------------------------------------------")
  
#-----------------------------------------------
# Timers for the code...

  timeUpdate  = myTimer(name = "Update")
  timePrint   = myTimer(name = "Print")
  timeIO      = myTimer(name = "I/O")
  timeMain    = myTimer(name = "MAIN Program", minutes=True)
  timeHxF     = myTimer(name = "HxF")
  
# Time the ESRF analysis

  timeMain.start()

# Parse input command options

  parser = OptionParser()
  parser.add_option("-o", "--obs",        dest="obs",       type="string",  help = "Observation file")
  parser.add_option("-t", "--time",       dest="time",      type="string",  help = "Time of prior calculations: 2008,05,8,22,10,00")
  parser.add_option("-e", "--exper",      dest="exper",     type="string",  help = "experiment run file created to store database info")
  parser.add_option(      "--nthreads",   dest="nthreads",  type="int",     help = "Number of threads for LETKF computation")
  parser.add_option(      "--window",     dest="window",    type="int",     help = "Width in seconds of Hx window")
  parser.add_option(      "--init",       dest="init",      default=False,  help = "Create new prior file", action="store_true")
  parser.add_option(      "--final",      dest="final",     default=False,  help = "Movie prior file to new DT-labeled file", action="store_true")
  parser.add_option(      "--obserr",     dest="obserr",    default=None,   type = "string", nargs=4, help = "Set observational errors, e.g., VR 3.0 DBZ 7.5")
  
  (options, args) = parser.parse_args()
  
#-------------------------------------------------------------------------------  
# Get experiment file

  if options.exper == None:
    print("\n --> ComputeHx:  No experiment's run filename specified, exiting....\n")
    parser.print_help()
    sys.exit(-1)
  else:
    with open(options.exper, 'rb') as p:
        exper = json.load(p)
        
# get path for file creation and location

  path = exper['base_path']
  
#-------------------------------------------------------------------------------
# Get the time stamp
    
  if options.time == None:
    print("\n --> ComputeHx:  No analysis time specified, exiting....\n")
    parser.print_help()
    sys.exit(-1)
  else:
    time = options.time.split(",")

#-------------------------------------------------------------------------------
# Obs error

  if options.obserr == None:
    obs_error = exper['DA_PARAMS']['obs_errors']
    print("\n --> ComputeHx:  using the obs errors from the EXPER file:  %s" % obs_error)
  else:
    print("\n --> ComputeHx:  using the obs errors from the command line")
    obs_error = _default_obs_error
    for key in obs_error.keys():
      for n, item in enumerate(options.obserr):
        if item == obs_error[key][0]:
           print("\n --> ComputeHx:  Changing %s observational error to:  %s" % (obs_error[key][0],options.obserr[n+1]))
           obs_error[key][1] = float(options.obserr[n+1])

#-------------------------------------------------------------------------------
# Set up window to look for

  if options.window == None:
    window = int(exper['DA_PARAMS']['assim_window'])
    print("\n --> ComputeHx:  using the window from the EXPER file:  +/- %d secs" % (window/2))
  else:
    window = int(options.window)
    print("\n --> ComputeHx:  using the window from the command line:  +/- %d secs" % (window/2))

#-------------------------------------------------------------------------------    
# Observation file

  if options.obs == None:
    obs_file = exper['radar_obs']
    print("\n --> ComputeHx:  using the observation file from the EXPER file:  %s" % obs_file)
  else:
    obs_file = options.obs
    print("\n --> ComputeHx:  using the observation file from the command line:  %s" % obs_file)

#-------------------------------------------------------------------------------    
# Init flag creates a new file

  if options.init:
    init = options.init
  else:
    init = False

#-------------------------------------------------------------------------------    
# Final flag moves prior file to Prior_DateTime.nc

  if options.final:
    final = options.final
  else:
    final = False

#-------------------------------------------------------------------------------    
# Need this value for 

  outlier = exper['DA_PARAMS']['outlier']
  print("\n --> ComputeHx: obs outlier is from EXPER file:  %d" % outlier)

   
#################################################################################################
# Read and search the observation file for each chunk of time...

  timeHxF.start()
  
  ob_f = pyDart.pyDART()
  ob_f.file(obs_file)   
  
# Initialize a bunch of containers, they get converted to numpy arrays below

  lat     = N.empty((0))
  lon     = N.empty((0))
  dates   = N.empty((0))
  value   = N.empty((0))
  kind    = N.empty((0))
  height  = N.empty((0))
  elev    = N.empty((0))
  az      = N.empty((0))
  err_var = N.empty((0))
  idx     = N.empty((0))

  analysis_time = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))
     
  print "\n --> ComputeHx:  Reading in model state for HxF at %s \n" %  (analysis_time.strftime("%Y-%m-%d %H:%M:%S"))
    
# read from history files

  timeIO.start()
    
  files = ens.FindRestartFiles(exper, analysis_time, ret_exp=False, ret_DT=False)

  state = ens.read_CM1_ens(files, exper, state_vector=state.Hxf, DateTime=analysis_time)   
  ens.ens_CM1_C2A(state)
  ens.ens_CM1_coords(state)

  timeIO.stop()
    
  Hxfm = N.empty((0,state.ne))  # need to init this here cause need size of ensemble
    
  g_lat_max = state.late[:].max()
  g_lat_min = state.late[:].min()
  g_lon_max = state.lone[:].max()
  g_lon_min = state.lone[:].min()
  g_alt_max = max(state.zc.data[:]) + state.hgt
  g_alt_min = min(state.zc.data[:]) + state.hgt
    
  calt = " & ( " + str(g_alt_min) + " < height < " + str(g_alt_max) + " ) "
  clat = "( " + str(g_lat_min) + " < lat < " + str(g_lat_max) + " ) "
  condition = clat + " & ( " + str(g_lon_min) + " < lon < " + str(g_lon_max) + " ) "

  dt     = DT.timedelta(0,int(window/2))
  begin  = analysis_time - dt
  ending = analysis_time + dt
    
  print "\n --> ComputeHx:  Using pyDart to search with condition: " + condition  
  print "\n --> ComputeHx:  Using pyDart to search with begin time of: ", begin.strftime("%Y-%m-%d %H:%M:%S")
  print "\n --> ComputeHx:  Using pyDart to search with end   time of: ", ending.strftime("%Y-%m-%d %H:%M:%S")
  print begin.timetuple()[:6]
  print ending.timetuple()[:6]
    
# ob_f.search(start=begin.timetuple()[:6], end=ending.timetuple()[:6], condition=condition)
  ob_f.search(start=begin.timetuple()[:6], end=ending.timetuple()[:6])

# number of observations, if there are none, kick out of the loop...
  
  if len(ob_f.index) > 0:
    print "\n --> ComputeHx:  Total number of obs found at search time: %s \n" % len(ob_f.index)
    print ob_f.index
  else:
    print "\n --> ComputeHx:  No obs found at search time:  %s exiting......\n" % (analysis_time)
    sys.exit(0)
  
# using the search index generated from above, obtain the location, data, and type

  subdata = ob_f.get_data()

# Compute Hxfs from 

  idx, Hxf, kind, lat, lon, height, elev, azimuth = ens.calcHx(state, 
                                                               subdata['kind'], 
                                                               subdata['lat'],                                         
                                                               subdata['lon'],
                                                               subdata['height'],
                                                               subdata['elevation'], subdata['azimuth'])
   
# At this point we have created all the Hxfs and so we can make sure we have enough obs to run

  if idx != None:

    nobs = N.size(idx) 

    print "\n --> ComputeHx:  Total number of obs found is: %d \n" % (nobs)

# retrieve these from pyTable

    err_var = N.append(err_var, subdata['error_var'][idx])   
    value   = N.append(value,   subdata['value'][idx])
    dates   = N.append(dates,   subdata['date'][idx])
  
  else:
    print "\n  >======================================================================<\n"    
    print "\n    ComputeHx: NOBS == 0:  NO HXF's were created - EXITING computeHx!!!  "
    print "\n    ComputeHx:  NOBS == 0:  NO HXF's were created - EXITING computeHx!!! "
    print "\n  >======================================================================<\n"   
    sys.exit(0)
    
# Overide the file observation std deviations with either defaults at top of script or input parameters

  for n, kk in enumerate(kind):
    if obs_error.has_key(kk):  err_var[n] = obs_error[kk][1]**2.0
  
# Here we choose to create a coordinate system of x/y's based from the SW corner (lat,lon) of model grid.
# The observations' new x/y's are then in the model's coordinate system relative to reference (lat,lon) of grid

  xs, ys = pyDart.dll_2_dxy(state.late[0], lat, state.lone[0], lon, degrees=True)
  zs     = height
  
  Hxfbar = N.average(Hxf,axis=1)
    
# We dont use this flag (is supposed to be a flag to tell one whether the obs is used for verification or assimilation, or thrown out)
  status = N.ones(nobs)
  
# WRITE OUT priors

# Here we create the Prior file if needed.  Note, the prior file may not be same datetime as current time
  if init:
    create_prior_file(state.ne)

  write_prior(state.ne, kind, value, dates, err_var, xs, ys, zs, Hxf, Hxfbar, lat, lon, elev, azimuth, status, outlier)
  
  if final:
    newPriorFile = os.path.join(path, "Prior_%s.nc" % analysis_time.strftime("%Y-%m-%d_%H:%M:%S"))
    os.rename("Prior.nc", newPriorFile)
  
  timeHxF.stop()
  
  timeIO.printit("--> ComputeHx:  Total time for ensemble I/O")
  
  timeHxF.printit("--> ComputeHx:  Time for searching ob table and computing priors")

  timeMain.stop()
  timeMain.printit("--> ComputeHx:  Total time in minutes")

# Print out Wallclock time for ComputeHx

  now = DT.datetime.now()

  print("\n  ----------------------------------------------------------------------")
  print("\n                END PROGRAM ComputeHx                                   ")
  print("\n      WALLCLOCK END TIME:  %s " % now.strftime("%Y-%m-%d %H:%M")     )
  print("\n  ----------------------------------------------------------------------\n")

