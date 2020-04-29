import sys
import os
import numpy as N
import datetime as DT 
from netCDF4 import *
import pyDart
from pyDart import dll_2_dxy
import ens 
import glob
from optparse import OptionParser
from multiprocessing import Pool
from time import time as timer
from nc_prior_file import *
from assim_util import *
import state_vector as state
import pickle

# Fortran functions
sys.path.append( "./fsrc" )
from fpython2 import fstate

#-------------------------------------------------------------------------------
# Default values 

_missing_value = 9900000000.
_time_units    = 'seconds since 1970-01-01 00:00:00'
_calendar      = 'standard'

debug  = False
_update_theta = 1

#-------------------------------------------------------------------------------
# Run parameters for the filter

obs_diag = {11: ["VR"], 12:["DBZ"]}


#-------------------------------------------------------------------------------
def getIndexEqual(field, value):
    return ( field == value )

#-------------------------------------------------------------------------------
def getIndexGT(field, value):
    return ( field > value )

#-------------------------------------------------------------------------------
def getIndexLT(field, value):
    return ( field < value )

########################################################################################################################
#
# Main program

if __name__ == "__main__":

  stime = timer()
  now = DT.datetime.now()

  print "----------------------------------------------------------------------\n"
  print "              BEGIN PROGRAM LETKF                                     \n "
  print "  WALLCLOCK START TIME:  %s \n" % now.strftime("%Y-%m-%d %H:%M")  
  print "  --------------------------------------------------------------------\n"
  
# Create timers for the code...

  timeFilter   = myTimer(name = "Filter")
  timeIO       = myTimer(name = "I/O")
  timeMain     = myTimer(name = "Program")
  timeHxF      = myTimer(name = "HxF")
  
# Time the data assimilation

  timeMain.start()

#-------------------------------------------------------------------------------
# Parse input command options

  parser = OptionParser()
  parser.add_option("-f", "--file",      dest="file",     type="string", help = "Priors netCDF file")
  parser.add_option("-t", "--time",      dest="time",     type="string", help = "Analysis time in YYYY,MM,DD,HH,MM,SS")
  parser.add_option("-e", "--exper",     dest="exper",    type="string", help = "experiment run file created to store database info")
  parser.add_option(      "--nthreads",  dest="nthreads", type="int",    help = "Number of threads for LETKF computation")
  parser.add_option(      "--outl",      dest="outl",     type="int",    help = "Outlier threshold for observations")
  parser.add_option(      "--saveW",     dest="saveW",    default=False, help = "Save LETKF weights...", action="store_true")
  parser.add_option(      "--readW",     dest="readW",    default=False, help = "Read LETKF weights...", action="store_true")
  parser.add_option(      "--noupdate",  dest="noupdate", default=False, help = "Do not update variables", action="store_true")
  parser.add_option(      "--aInflate",  dest="ainflate",  type = "int",   help = "Sets type of adaptive inflation (0,1,2,3)", default=None)
  
  (options, args) = parser.parse_args()

#-------------------------------------------------------------------------------  
# Get experiment file

  if options.exper == None:
    print("\n --> LETKF:  No experiment's run filename specified, exiting....\n")
    parser.print_help()
    sys.exit(-1)
  else:
    with open(options.exper, 'rb') as p:
        exper = pickle.load(p)

#-------------------------------------------------------------------------------
# Get the time stamp
    
  if options.time == None:
    print("\n --> LETKF:  No analysis time specified, exiting....\n")
    parser.print_help()
    sys.exit(-1)
  else:
    time = options.time.split(",")
    
#-------------------------------------------------------------------------------
# The rest of the input is options - so set parameters based on what is in the experiment data base

  mpass             = exper['DA_PARAMS']['mpass']
  writeFcstMean     = exper['DA_PARAMS']['writeFcstMean']
  writeAnalMean     = exper['DA_PARAMS']['writeAnalMean']
  saveWeights       = exper['DA_PARAMS']['saveWeights']
  readWeights       = exper['DA_PARAMS']['readWeights']
  rhoriz            = exper['DA_PARAMS']['rhoriz']
  rvert             = exper['DA_PARAMS']['rvert']
  rtime             = exper['DA_PARAMS']['rtime']
  cutoff            = exper['DA_PARAMS']['cutoff']
  zcutoff           = exper['DA_PARAMS']['zcutoff']
  inflate           = exper['DA_PARAMS']['inflate']
  print_state_stats = exper['DA_PARAMS']['print_state_stats']
  path              = exper['base_path']

#-------------------------------------------------------------------------------
# Parallel Threading

  if options.nthreads != None:
    nthreads = options.nthreads
    print("\n --> LETKF:  Number of threads requested is from command line:  %d" % nthreads)
  else:
    nthreads = exper['DA_PARAMS']['nthreads']
    print("\n --> LETKF:  Number of threads requested is from EXPER file:  %d" % nthreads)
    

#-------------------------------------------------------------------------------
#  Priors file
    
  if options.file == None:
    print("\n --> LETKF:  No netCDF Priors file specified, using default priors file")
    ftime      = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))
    prior_file = os.path.join(path, "Prior_%s.nc" % ftime.strftime("%Y-%m-%d_%H:%M:%S"))
    if os.path.isfile(prior_file):
      print("\n --> LETKF:  Using prior file:  %s\n" % prior_file)
    else:
      print("\n  --> LETKF:  Could not find prior file:  %s !!!!  Exiting LETKF...." % prior_file)
      parser.print_help()
      sys.exit(-1)
  else:
    prior_file = options.file

#-------------------------------------------------------------------------------
# Inflation type
  if options.ainflate:
    aInflate = options.ainflate
  else:
    aInflate   = exper['DA_PARAMS']['aInflate']
    
  if aInflate == 0:
    inflate_file_exists = False
  else:
    inflation_deltaT = exper['DA_PARAMS']['assim_window']
    dt = DT.timedelta(0,inflation_deltaT)
    in_time = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5])) - dt
    inflate_file = os.path.join(path, "Inflation_%s.nc" % in_time.strftime("%Y-%m-%d_%H:%M:%S"))
    print("\n --> LETKF:  Trying to find inflation file: %s" % inflate_file)
    if os.path.isfile(inflate_file):
      inflate_file_exists = True
      print("\n --> LETKF:  Using inflation file: %s" % inflate_file)
    else:
      inflate_file_exists = False
      print("\n --> LETKF:  Could not find inflation file....none will be used.")

#-------------------------------------------------------------------------------
# LETKF Weight I/O flags

  if options.saveW:
    saveWeights = options.saveW
    print("\n --> LETKF: saveWeights is from command line:  %s" % saveWeights)
  else:  
    saveWeights = exper['DA_PARAMS']['saveWeights']
    print("\n --> LETKF: saveWeights is from EXPER file:  %s" % saveWeights)

  if options.readW:
    readWeights = options.readW
    print("\n --> LETKF: readWeights is from command line:  %s" % readWeights)
  else:  
    readWeights = exper['DA_PARAMS']['readWeights']
    print("\n --> LETKF: saveWeights is from EXPER file:  %s" % readWeights)
   
#-------------------------------------------------------------------------------
# Update flag which is useful when debugging (does not overwrite files

  if options.noupdate:
      noUpdate = options.noupdate
      print("\n --> LETKF: noUpdate is set from command line:  No files will be written back")
  else:
      noUpdate = False

#-------------------------------------------------------------------------------
# outlier threshold

  if options.outl == None:
    outlier_threshold = exper['DA_PARAMS']['outlier']
    print("\n --> LETKF: obs outlier is from EXPER file:  %d" % outlier_threshold)
  else:
    outlier_threshold = options.outl
    print("\n --> LETKF: obs outlier is from command line:  %d" % outlier_threshold)

########################################################################################################################
#
# Read in STATE VECTOR for ensemble
  
  timeIO.start()

  analysis_time = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))
  
  files, exper = ens.FindRestartFiles(exper, analysis_time, ret_DT=False)
  
  print "\n  ---->LETKF:  Reading in model state for analysis at %s \n" %  (analysis_time.strftime("%Y %m-%d %H:%M:%S"))
  
  state = ens.read_CM1_ens(files, exper, state_vector = None, DateTime=analysis_time, addmean=1)   
  ens.ens_CM1_coords(state)
  ens.ens_CM1_mean(state)
  
  tanalysis = date2num(analysis_time,units=_time_units,calendar=_calendar)
  
  print "\n  ---->Analysis state read in for time %s " % (analysis_time.strftime("%Y %m-%d %H:%M:%S"))
  print "\n  ---->Analysis state universal time %f" % (tanalysis)

  if debug:
    print("WARNING WARNING WARNING:  LETKF is in DEBUG mode!!!")
    print("WARNING WARNING WARNING:  LETKF is in DEBUG mode!!!")
    print("WARNING WARNING WARNING:  LETKF is in DEBUG mode!!!")
    print("WARNING WARNING WARNING:  No files are written!!!!")
  else:
    if writeFcstMean:
      ens.write_CM1_ens(state, writeFcstMean=writeFcstMean, overwrite=True)
  
  timeIO.stop()
  
########################################################################################################################
#
# Read in observational priors

  timeHxF.start()
  
  analysis_sec = date2num(analysis_time,units=_time_units,calendar=_calendar)
     
  print("\n  ----> LETKF:  Reading in Prior file:  %s \n" %  (prior_file))
    
  kind, value, error, lat, lon, xob, yob, zob, tob, Hxf, Hxfbar, outlier = read_prior(analysis_time, file=prior_file, sort=True)

  nobs = Hxf.shape[0]
  ne   = Hxf.shape[1]
  
  if nobs == 0:  
    print "\n  ----> LETKF:  Prior file is empty, exiting at analysis time %s \n" %  (analysis_time.strftime("%Y %m-%d %H:%M:%S"))
    sys.exit(0)
  else:
    print "\n  ----> LETKF:  Prior file has %d obs at analysis time %s \n" %  (nobs, analysis_time.strftime("%Y %m-%d %H:%M:%S"))

    tob = N.array(tob - analysis_sec, dtype=N.float64)   # departure time for temporation localization

    if (outlier_threshold > 0):
      print('\n  --> LETKF called with an outlier threshold of %d standard deviations' % outlier_threshold)

      mask  = N.where( outlier <= outlier_threshold, True, False )
      mask2 = N.where( kind == 11, True, False )
      print "Mask before Vr mask:  ", N.count_nonzero(mask)
      mask  = mask | mask2
      print "Mask after Vr mask:  ", N.count_nonzero(mask)

      kind   = kind[mask]
      value  = value[mask]
      error  = error[mask]
      tob    = tob[mask]
      xob    = xob[mask]
      yob    = yob[mask]
      zob    = zob[mask]
      Hxf    = Hxf[mask][:]
      Hxfbar = Hxfbar[mask]
      lat    = lat[mask]
      lon    = lon[mask]

      print('\n  >=====================================================================<')
      print('\n  --> LETKF:  Number of total    obs = %d' % nobs)
      print('\n  --> LETKF:  Number of rejected obs = %d' % (nobs - N.count_nonzero(mask) ))
      print('\n  --> LETKF:  Max outlier: %f  Min outlier:  %f\n' % (outlier.max(), outlier.min()))
    
# Reset number of obs

      nobs = N.count_nonzero(mask)

# Recompute coordinate system of x/y's based from the SW corner (lat,lon) of model grid.
# The observations' new x/y's are then in the model's coordinate system relative to reference (lat,lon) of grid
# Need to do this again because of running-in-place for moving grids...

  xob, yob = pyDart.dll_2_dxy(state.late[0], lat, state.lone[0], lon, degrees=True)
  dep      = N.zeros(nobs)
  dep[:]   = value[:] - Hxfbar[:]
    
# Here we create the localization arrays

  rdiag    = N.array(error.reshape(nobs,1), order='F')
  rloc     = N.array(N.ones(nobs).reshape(nobs,1), order='F')

# Print some stats...

  print "\n  >==============================================================================================<"
  print "\n  LETKF:  %s  Total number of observations:  %d" % (analysis_time.strftime("%Y-%m-%d_%H:%M:%S"),value.size)
  print "\n  <==============================================================================================<"
  print "\n  LETKF:  Obs Space Diagnostic for forecast computed via Dowell and Wicker, MWR 2009, Additive Noise\n"
  print "\n  LETKF:  Obs Std  |  Root-Mean-Sq-Error  |    Bias    |  Spread (obsErr^2+Hxf_var^2)  |   Consistency Ratio"
  print "\n  >==============================================================================================<\n"

  for key in obs_diag.keys():
    index_kind  = getIndexEqual(kind, key)
    d           = dep[index_kind]
    
    if obs_diag[key][0] == 'DBZ':

  # Non-zero dBZ obs
  
        index_dbz   = index_kind & getIndexGT(value, 0.1)
        if N.sum(index_dbz) > 0:
          ob_err      = error[index_dbz].mean()
          d           = dep[index_dbz]
          Hxftmp      = Hxf[index_dbz,:]
          Hxf_var     = Hxftmp.var(ddof=1, axis=1).mean()
          inno_var    = N.mean((d - d.mean())**2)
          consi_ratio = (ob_err + Hxf_var) / inno_var

          print("\n -->  LETKF:  %s  NOBS: %5.5d  %3.3s>0: %3.1f  RMSI: %6.3f  M-Innov: %7.3f  Spread: %6.3f  CRatio: %7.4f " \
                % (analysis_time.strftime("%Y-%m-%d_%H:%M:%S"), N.sum(index_dbz), obs_diag[key][0], N.sqrt(ob_err), \
                   N.sqrt(inno_var), d.mean(), N.sqrt(ob_err + Hxf_var), consi_ratio))
        else:
          print("\n -->  LETKF:  %3.3s>0: is not present in observations" % (obs_diag[key][0]))
        
  # zero dBZ obs      

        index_dbz   = index_kind & getIndexLT(value, 0.1)
        if N.sum(index_dbz) > 0:
          ob_err      = error[index_dbz].mean()
          d           = dep[index_dbz]
          Hxftmp      = Hxf[index_dbz,:]
          Hxf_var     = Hxftmp.var(ddof=1, axis=1).mean()
          inno_var    = N.mean((d - d.mean())**2)
          consi_ratio = (ob_err + Hxf_var) / inno_var
          print("\n --> LETKF:%s  NOBS: %5.5d  %3.3s<0: %3.1f  RMSI: %6.3f  M-Innov: %7.3f  Spread: %6.3f  CRatio: %7.4f " \
                % (analysis_time.strftime("%Y-%m-%d_%H:%M:%S"), N.sum(index_dbz), obs_diag[key][0], N.sqrt(ob_err), \
                   N.sqrt(inno_var), d.mean(), N.sqrt(ob_err + Hxf_var), consi_ratio))
        else:
          print("\n -->  LETKF:  %3.3s<=0: is not present in observations" % (obs_diag[key][0]))

  # other obs

    elif d.size > 0:

      ob_err      = error[index_kind].mean()
      d           = dep[index_kind]
      Hxftmp      = Hxf[index_kind,:]
      Hxf_var     = Hxftmp.var(ddof=1, axis=1).mean()
      inno_var    = N.mean((d - d.mean())**2)
      consi_ratio = (ob_err + Hxf_var) / inno_var
      print("\n -->  LETKF:  %s  NOBS: %5.5d    %3.3s: %3.1f  RMSI: %6.3f  M-Innov: %7.3f  Spread: %6.3f  CRatio: %7.4f " \
                % (analysis_time.strftime("%Y-%m-%d_%H:%M:%S"), N.sum(index_kind), obs_diag[key][0], N.sqrt(ob_err), \
                   N.sqrt(inno_var), d.mean(), N.sqrt(ob_err + Hxf_var), consi_ratio))
    else:
      print "\n -->  LETKF:  %s is not present in observations" % (obs_diag[key][0])

  print "\n  >============================================================================================<\n"

  timeHxF.stop()
  
#################################################################################################
#
# Adaptive inflation: read in or create dummy file....

  if aInflate == 0:   # just create a dummy array of ones...
    print("\n -->  LETKF:  Adaptive inflation is turned off....\n")
    inflate = N.ones((state.nz, state.ny, state.nx), order='F', dtype=N.float64)
  else:
    print("\n -->  LETKF:  Adaptive inflation is turned on....\n")
    if inflate_file_exists:
      print("\n -->  LETKF:  Adaptive inflation is using prior file:  %s\n" % inflate_file)
      inflate = read_inflation_file(file=inflate_file)
    else:
      print("\n -->  LETKF: Adaptive inflation cannot find the previous inflation file, initializing inflation array to 1.0 \n") 
      inflate = N.ones((state.nz, state.ny, state.nx), order='F', dtype=N.float64)
      
################################################################################################################################
#
# Begin LETKF analysis
 
  if print_state_stats: print_state_diagnostics(state, store=True)
  
  print("\n  >=======================  >>BEGIN FILTER<< ===========================================<\n")    

# Call the fortran letkf module

  timeFilter.start()

  if mpass == True:

    print("\n -->  LETKF:  MultiPASS algorithm turned on, assimilating VR THEN DBZ\n")
    for key in obs_diag.keys():

      index_kind  = getIndexEqual(kind, key)
      if obs_diag[key][0] == "DBZ":  
        update_theta = 0
      else:
        update_theta = _update_theta

      if N.sum(index_kind) > 0:
        print("\n-->  LETKF:  %d observations of %3.3s is now being assimilated!\n" % (N.sum(index_kind), obs_diag[key][0]))
        inflateN = fstate.compute_letkf(xob[index_kind], yob[index_kind], zob[index_kind], tob[index_kind], tanalysis,
                                        value[index_kind], Hxf[index_kind], dep[index_kind], rdiag[index_kind], rloc[index_kind], 
                                        rhoriz, rvert, rtime, nthreads, cutoff, zcutoff, update_theta, aInflate, inflate,
                                        saveWeights, readWeights)
      else:
        print("\n -->  LETKF:%3.3s: NO OBSERVATIONS FOUND!\n" % (obs_diag[key][0]))

  else:

    print("\n -->  LETKF:ALL OBS TYPES BEING ASSIMULATED SIMULTANEOUSLY: MULTIPASS==False\n")
    update_theta = _update_theta
    inflateN = fstate.compute_letkf(xob, yob, zob, tob, tanalysis,
                                    value, Hxf, dep, rdiag, rloc, 
                                    rhoriz, rvert, rtime, nthreads, cutoff, zcutoff, update_theta, aInflate, inflate,
                                    saveWeights, readWeights)
                           
# inflateN always has the adaptive inflation field produced by LETKF core. Write it out...

  if aInflate > 0:
    print("\n -->  LETKF:  Adaptive Inflation stats:  Max AI:  %3.1f  Min AI %3.1f  Mean AI: %3.1f  StdDev:  %4.3f\n" % \
        (inflateN.max(), inflateN.min(), inflateN.mean(), inflateN.std()))

  write_inflation_file(inflateN,fstate.xc,fstate.yc,fstate.zc,analysis_time)

# If saveWeights == True, rename the netCDF4 file dumped out by the fortran code

  if saveWeights:
    rename_weight_file(analysis_time)

  timeFilter.stop()

  print "\n  >=======================  >>END FILTER<< ===========================================<\n"    

  if print_state_stats: print_state_diagnostics(state, header=True)   
  
# Write out model state
  
  timeIO.start()

# Write ensemble back out

  if debug:
    print("\nWARNING WARNING WARNING:  LETKF is in DEBUG mode!!!")
    print("\nWARNING WARNING WARNING:  LETKF is in DEBUG mode!!!")
    print("\nWARNING WARNING WARNING:  LETKF is in DEBUG mode!!!")
    print("\nWARNING WARNING WARNING:  STATE DATA IS NOT WRITTEN BACK OUT!")
  elif noUpdate:
    print("\n -->  LETKF:  NO UPDATE FLAG is TRUE....state data will not be written back out!")
  else:
    print("\n -->  LETKF:  state data now written back out!")
    ens.write_CM1_ens(state, writeEns=True, writeAnalMean=writeAnalMean, overwrite=False)
  
  timeIO.stop()
  
# Print out specific timing information

  timeHxF.printit("")
  timeFilter.printit("")
  timeIO.printit("")
  
  timeMain.stop()
  timeMain.printit("")

# Print out Wallclock time for LETKF

  now = DT.datetime.now()

  print "\n  ----------------------------------------------------------------------"
  print "\n                END PROGRAM LETKF                                        "
  print "\n  WALLCLOCK END TIME:  %s \n" % now.strftime("%Y-%m-%d %H:%M")  
  print "\n  -----------------------------------------------------------------------\n"
