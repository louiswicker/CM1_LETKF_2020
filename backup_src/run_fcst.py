#!/usr/bin/env python
#
# System imports
#   
import os, sys, glob
import time as cpu
import string
from optparse import OptionParser
import numpy as N
from multiprocessing import Pool
import f90nml
import pickle
import netCDF4 as ncdf
import datetime as DT

_nthreads = 20
_restart_frequency = 300.

debug = False

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

#=======================================================================================================================
# RunMember is a function that runs a system command

def RunMember(cmd):
    print("\n Executing command:  %s " % cmd)
    os.system(cmd)
    print("\n %s completed...." % cmd)
    return
    
#-------------------------------------------------------------------------------
#
# Command line arguments

usage  = "usage: %prog [options] arg"

parser = OptionParser(usage)
parser.add_option("-e", "--exp",      dest="exp",     type="string", help = "Path to the pickled experiment file generated \
                                                                                 from the create_run_letkf script")
parser.add_option("-t", "--time",     dest="datetime", default=None, type = "string", help = "Usage:  --time 2003,5,8,21,0,0")   

parser.add_option(      "--run_time", dest="run_time", type="int",    help = "Run time of model forecast in seconds")
parser.add_option(      "--nthreads", dest="nthreads", type="int",    help = "Number of threads to run model")
parser.add_option("-i", "--init",     dest="init",     default=False, action="store_true", help = "Create initial conditionj file")
parser.add_option(      "--range",    dest="range",    type ="int",   default=None, nargs=2, help = "The range of ensemble members to run forecasts \
                                                                                 for. The default is all members. Usage:  --range 1 36")

(options, args) = parser.parse_args()

if options.exp == None:
    parser.print_help()
    print "\n RUN_FCST: ERROR -->  Experiment file not supplied..EXITING!!!"
    sys.exit(0)
else:
    with open(options.exp, 'rb') as f:
        experiment = pickle.load(f)

if options.init == True:
    init = True
else:
    init = False
    if options.datetime == None:
         parser.print_help()
         print("\n ==> RUN_FCST: ERROR --> date and time for start time not supplied..EXITING!!!")
         sys.exit(-1)
    else:
         list = [int(t) for t in options.datetime.split(",")]
         startDT    = DT.datetime(list[0],list[1],list[2],list[3],list[4],list[5])
         runDT      = DT.datetime(experiment["YEAR"],experiment["MONTH"],experiment["DAY"],experiment["HOUR"],experiment["MINUTE"],experiment["SECOND"])
         start_time = (startDT-runDT).seconds
         print("\n ==> RUN_FCST: Date&time supplied: %s, which corresponds to a local model time of %d sec" % (startDT.strftime("%Y %m-%d %H:%M:%S"), start_time))

if options.run_time == None:
    print "\n Run_Fcst Script: IMPORTANT:  Run time not supplied...defaulting to init mode!!!"
    options.init = True
    run_time = 0.0
else:
    run_time = options.run_time

if options.nthreads != None:
    nthreads = options.nthreads
    print("\n Run_Fcst Script:  %d threads requested...." % nthreads)
else:
    nthreads = _nthreads
    print("\n Run_Fcst Script:  defaulting to %d threads" % nthreads)
    
if options.range == None:
    ne = experiment['ne']
    ne_start = 1
    ne_end   = ne
    print("\n Run_Fcst Script: Number ensemble members to run not supplied, defaulting to experiment file:  %d" % ne)
else:
    ne = 1 + options.range[1] - options.range[0]
    ne_start = options.range[0]
    ne_end   = options.range[1]
    print("\n Run_Fcst Script: Number ensemble members supplied  %d  %d" % ne_start, ne_end)

#-------------------------------------------------------------------------------
#
# Integrate ensemble members to next observation time.
#
cpu_model = 0.0
c0 = cpu.time()

if init:
    print("\n Run_Fcst Script: Initializing %3.3d ensemble members" % (ne))
else:
    print("\n Run_Fcst Script: Integrating %3.3d ensemble members for %4.4d seconds" % (ne, run_time))

nthreads = min(nthreads, ne)

#-----------------------------------------------------------------------------------------------------
# Write now we assume that the namelist each member uses is the same namelist.  So we just need to 
# edit the namelist.input file found at the top level of the run directory.  
# 
# HANDY Python module:  f90nml  -->  this is why python is great, someone already did this!!!!!!

namelist = f90nml.read(os.path.join(experiment['base_path'],"namelist.input"))

if init:
    namelist['param1']['run_time']         = 0
    namelist['param2']['irst']             = 0
    namelist['param1']['rstfrq']           = 0.0
    namelist['param9']['restart_format']   = 2
    namelist['param9']['restart_filetype'] = 1
    namelist['param2']['rstnum']           = 0

else:
    output_basename = namelist['param9']['output_basename']

    namelist['param2']['irst']             = 1
    namelist['param2']['rstnum']           = FindRestartFile(experiment['fcst_members'][0],output_basename,start_time)
    namelist['param1']['run_time']         = run_time
    namelist['param9']['restart_format']   = 2
    namelist['param9']['restart_filetype'] = 1

# assim_freq > 0 tells you that you need model restart files dumped every freq secs so that you   
# can compute the priors closer to the observations for 4D LETKF - e.g.. enables asynchronous DA

    if experiment['DA_PARAMS']['assim_freq'] > 0:   
        namelist['param1']['rstfrq'] = experiment['DA_PARAMS']['assim_freq']
        
# assim_freq < 0 tells you that the data assimilation is synchronous
# Only need file model restart files at the frequency of the assimilation window.

    else:                                           
        namelist['param1']['rstfrq'] = experiment['DA_PARAMS']['assim_window']
    
# Write out new namelist for runs

namelist.write(os.path.join(experiment['base_path'],"namelist.input"), force=True)

#-----------------------------------------------------------------------------------------------------
#
# Now actually run the forecasts
#

pool = Pool(processes=nthreads)              # set up a queue to run

for n in N.arange(ne_start-1,ne_end):

    fcst_member = experiment['fcst_members'][n]
    model       = os.path.join(fcst_member, "cm1.exe")
    outputfile  = os.path.join(fcst_member, "cm1.out")

    if debug:  
        print("\n Now setting up forecast member: %s" % fcst_member)
        print("%s is the model path" % model)
        print("%s is the outputfile path" % outputfile)

    cmd = "cd %s ; %s >& %s" % (fcst_member, "cm1.exe", "cm1.out")
    pool.apply_async(RunMember, (cmd,))

pool.close()
pool.join()

cpu_model = cpu_model + cpu.time() - c0

if init:
    print "\nRun_Fcst Script: Model runs initialized for experiment   %f  CPU secs\n" % (round(cpu_model, 3))
else:
    print "\nRun_Fcst Script: Model runs integrated out %d seconds, used   %f  CPU secs\n" % (run_time, round(cpu_model, 3))
