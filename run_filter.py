#!/usr/bin/env python
import sys
import os
import numpy as N
import datetime as DT 
from optparse import OptionParser
import glob
import subprocess
import pickle 
import json
########################################################################################################################
#
# This helps return the correct pipes from the parallel threading...

def run_unix_cmd(cmd, file_handle=None):
    if file_handle == None:
      p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
      return p.communicate()[0]
    else:
      p = subprocess.Popen(cmd, shell=True, stdout=file_handle, stderr=subprocess.STDOUT)
      return p.communicate()[0]

########################################################################################################################
#
# Main program

if __name__ == "__main__":

  print("\n  ----------------------------------------------------------------------")
  print("\n               BEGIN PROGRAM RUN_FILTER                                 ")
  print("\n  This script sets the parameters to create Hx's and runs the letkf code")
  print("\n  ----------------------------------------------------------------------\n")
  
# Parse input command options

  parser = OptionParser()
  parser.add_option("-o", "--obs",        dest="obs",       type="string",  help = "radar observation file in H5 (pyDart) format")
  parser.add_option("-t", "--time",       dest="time",      type="string",  help = "Analysis time in YYYY,MM,DD,HH,MM,SS")
  parser.add_option("-e", "--exper",      dest="exper",     type="string",  help = "experiment run file created to store database info")
  parser.add_option(      "--window",     dest="window",    type="int",     help = "Window (in sec) for assimilation, assumed +/- window/2")
  parser.add_option(      "--freq",       dest="freq",      type="int",     help = "Frequency (in secs) of history files in window, \
                                                                                    if freq=width, then DA is synchronous")
  parser.add_option(    "--nthreads",     dest="nthreads",  type="int",     help = "Number of threads for LETKF computation")
  parser.add_option(    "--aInflate",     dest="ainflate",  type = "int",   help = "Sets type of adaptive inflation (0,1,2,3)", default=None)
  parser.add_option(    "--HF",           dest="HF",        type = "int",   help = "Sets up writeback skip for HF filter", default=None)
  parser.add_option(      "--obserr",     dest="obserr",    type="string",  nargs=4, help = "")

  (options, args) = parser.parse_args()

#-------------------------------------------------------------------------------  
# Get experiment file

  if options.exper == None:
    print("\n --> Run_Filter:  No experiment's run filename specified, exiting....\n")
    parser.print_help()
    sys.exit(-1)
  else:
    with open(options.exper, 'rb') as p:
        #        exper = pickle.load(p)
        exper = json.load(p)

# get path for file creation and location

  path = exper['base_path']
  
#-------------------------------------------------------------------------------
# Get the time stamp
    
  if options.time == None:
    print("\n --> Run_Filter:  No analysis time specified, exiting....\n")
    parser.print_help()
    sys.exit(-1)
  else:
    time = options.time.split(",")

#-------------------------------------------------------------------------------
# Set up data assimilation window

  if options.window == None:
    window = int(exper['DA_PARAMS']['assim_window'])
    print(("\n --> Run_Filter:  using the window from the EXPER file:  +/- %d secs" % (window/2)))
  else:
    window = int(options.window)
    print(("\n --> Run_Filter:  using the window from the command line:  +/- %d secs" % (window/2)))

#-------------------------------------------------------------------------------
# Set up frequency of assimilation (needed for asynchronous mostly)

  if options.freq == None:
    freq = int(exper['DA_PARAMS']['assim_freq'])
    if freq > 0:
      print(("\n --> Run_Filter:  using Asynchronous DA with freq from EXPER file:  %d secs" % (freq)))
    else:
      print(("\n --> Run_Filter:  using Synchronous DA with freq from EXPER file:  %d secs" % (freq)))
  else:
    freq = options.freq
    if freq > 0:
      print(("\n --> Run_Filter:  using Asynchronous DA with freq from command line:  %d secs" % (freq)))
    else:
      print(("\n --> Run_Filter:  using Synchronous DA with freq from command line:  %d secs" % (freq)))

#-------------------------------------------------------------------------------
# Observation file

  if options.obs == None:
    obs_file = exper['radar_obs']
    print(("\n --> Run_Filter:  using the observation file from the EXPER file:  %s" % obs_file))
  else:
    obs_file = options.obs
    print(("\n --> Run_Filter:  using the observation file from the command line:  %s" % obs_file))

#-------------------------------------------------------------------------------
# Parallel Threading

  if options.nthreads != None:
    nthreads = options.nthreads
    print(("\n --> RUN_Filter:  Number of threads requested is from command line:  %d" % nthreads))
  else:
    nthreads = exper['DA_PARAMS']['nthreads']
    print(("\n --> Run_Filter:  Number of threads requested is from EXPER file:  %d" % nthreads))

#-------------------------------------------------------------------------------
# Inflation type

  if options.ainflate:
    aInflate = options.ainflate
    print(("\n --> Run_Filter:  using the inflation method from command line:  %d" % (aInflate)))
  else:
    aInflate = exper['DA_PARAMS']['aInflate']
    print(("\n --> Run_Filter:  using the inflation method from EXPER file:  %d" % (aInflate)))

#################################################################################################
#
# CREATE Hx's.  
#
# Assumption here is that one can sub-divide the input into N-frequency minute bins cleanly, 
#            this assumes 1, 3, 5, 7, ... minute windows
#    
  if freq > 0:
    h_width = ((int(window/freq))/2)*freq

# If freq is less than zero, that is a flag to use only one history file for Hxfs.  Essentially, the code runs synchronous
  else:
     h_width = 0
     freq    = -freq
     
  print("\n  >=======================  Run_Filter:  Now creating Priors =============================<  \n")
  
  for n, sec in enumerate(N.arange(-h_width,h_width+freq,freq)):
  
    dt      = DT.timedelta(seconds=N.int(sec))
    Hx_time = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5])) + dt
    print("TESTE")
    print(Hx_time)
    print(("\n  -->  Run_Filter calling computeHx for time %s  \n" % Hx_time.strftime("%Y,%m,%d,%H,%M,%S")))

    cmd = "computeHx.py --exper %s --time %s -o %s --window %d" % (options.exper, Hx_time.strftime("%Y,%m,%d,%H,%M,%S"), obs_file, freq)
                               
    if n == 0:  cmd = "%s --init" % cmd

    if options.obserr:  cmd = "%s --obserr %s %s %s %s" % (cmd, \
                        options.obserr[0],options.obserr[1],options.obserr[2],options.obserr[3])
       
    print(("\n  "+cmd+"\n"))

    print((run_unix_cmd(cmd)))

  file_DT = DT.datetime(int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))

  newPriorFile = os.path.join(path, "Prior_%s.nc" % file_DT.strftime("%Y-%m-%d_%H:%M:%S"))

  os.rename("Prior.nc", newPriorFile)
  
  print(("\n  --> Moved Prior.nc file to %s\n" % (newPriorFile)))

#################################################################################################
# Run LETKF 
#

  print(("\n  >================ Run_Filter:  Running LETKF at time: %s ===================<  \n" % file_DT.strftime("%Y-%m-%d_%H:%M:%S")))

  cmd = "python letkf.py --exper %s --time %s --nthreads %d " % (options.exper, options.time, nthreads)
  
  if aInflate:   cmd = "%s --aInflate %d" % (cmd, aInflate)
  if options.HF != None: cmd = "%s --HF %d" % (cmd, options.HF)
  
  print(("\n  "+cmd+"\n"))
  print((run_unix_cmd(cmd)))

  print("\n  ----------------------------------------------------------------------")
  print("\n                 END PROGRAM RUN_FILTER                                 ")
  print("\n  ----------------------------------------------------------------------")
