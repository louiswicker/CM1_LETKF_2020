#!/usr/bin/env python
import sys
import glob
import os
import time as cpu
import pickle

import numpy as N
from optparse import OptionParser
import datetime as DT
from subprocess import *
import json

#-----------------------------------------------------------------------------------------------------------------------------
# Define initialization strings - not needed if you run them yourself

ic_cmds = ["create_run_letkf.py",
           "run_fcst.py -e RUN_LETKF/RUN_LETKF.exp -i",
           "ens.py -e RUN_LETKF/RUN_LETKF.exp --init0 -t 2003,5,8,20,40,0 --write"]

#-----------------------------------------------------------------------------------------------------------------------------
#
def run_unix_cmd(cmd):
        p = Popen(cmd, shell=True, stdout=PIPE,stderr=PIPE)
        output, errors = p.communicate()
        return output, errors
#
#-------------------------------------------------------------------------------
#
# Command line arguments

usage  = "usage: %prog [options] arg"
parser = OptionParser(usage)

parser.add_option("-e",  "--exp", dest="exp", type="string", help = "Path to the pickled experiment file generated \
                                                                     from the create_run_letkf script")
                                                                     
parser.add_option("-d",  "--dir", dest="dir", type="string", help = "Experiment directory - assumes that the experiment \
                                                                     database file is named the same as directory")

parser.add_option("-i",  "--init", dest="init", default=False, action="store_true", help = "Run the initial setup commands specified \
                                                                     at the top of the file")

parser.add_option(  "--debug", dest="debug", default=False, action="store_true", help = "Setting this flag will dump the commands\
                   that would be run by the script - and actually, the comands could be dumping into a shell script\
                   for editing and running")

(options, args) = parser.parse_args()

if options.debug == True:
    RunIt = False
else:
    RunIt = True      # False just prints out what it will do, good for debugging
#
#-----------------------------------------------------------------------------------------------------------------------------
# Run anly defined initialization strings - not needed if you run them yourself

if options.init:
    cwd = os.getcwd()
    for cmd in ic_cmds:
        newcmd = os.path.join(cwd, cmd)
        if RunIt:
            print(("Running  %s" % newcmd))
            os.system(newcmd)
        #   print(("\n ==> run_Exper: ic_cmds:  STDOUT: %s " % master_output))
        #   print(("\n ==> run_Exper: ic_cmds:  STDERR: %s " % master_error))
        else:
           print(("%s" % newcmd))

#----------------------------------------------------------------------------------------------------------------------------
# Observation files:  (the ObFile can have more than one data file, add to list)

ObFile   = ["Obs/obs_seq_PAR_4km_1min.h5"]      # HDF5 file for observations to be used for assimilation (not used currently)

newstart = [False, 2003,5,8,20,40,0 ]
stop     = [2003,5,8,21,30,0]

fcst     = 900 

HF       = [False,2]
assim_count = 0

#----------------------------------------------------------------------------------------------------------------------------
# da params is simply a subset of what is in the create-letkf-run....
#
# NOT USED TO OVERWRITE CREATE_LETKF_RUN file as of YET!!!!!
# THIS SECTION DOES NOTHING, but code could be added below to "overwrite" the settings in the experiment database file.

da_params = {
             "obs_errors":     {           # the errors that are used for observations
                                11:  ["VR", 3.0],
                                12:  ["DBZ",7.5],
                               },
             "aInflate":            1,        # type(int): 0 => no adapt inflat / 1 => LETKF AI / 2 => WH2010 RTPS / 3 => RTPP
             "outlier":             3,        # type(int): Outlier threshold:  None means dont threshold, else set to sigma (e.g., 3, 5, 7 etc.)
             "nthreads":           12,        # type(int):  number of threads used to run the ensemble members and enkf (if parallel)
             "assim_window":      300,        # type(int):  window for assimilation
                                              #             (note that the assim window will be +/- (assim_window/2) )
             "assim_freq":       -300,        # type(int):  used to set asynchronous DA assimilation
             "cook":              600,        # type(int): time to pre-cook initial perturbations
            }

#-------------------------------------------------------------------------------
# 

if options.dir and options.exp == None:
    options.exp = glob.glob(os.path.join(options.dir, "*.exp" ))[0]
    print(("\n#==> run_Exper: found experiment files %s" % options.exp))

if options.exp == None:
    parser.print_help()
    print("\n ==> run_Exper: ERROR --> Experiment file not supplied..EXITING!!!")
    sys.exit(-1)
else:
    with open(options.exp, 'rb') as f:
        exper = json.load(f)


cpu_enkf     = 0.0
cpu_model    = 0.0
cpu_correct  = 0.0
cpu_total    = 0.0

#-------------------------------------------------------------------------------
#
# Create blocks of time needed to run all the assimilation....
# Create a master list of DT objects that are used to cycle the system.
#
assim_freq  = exper['DA_PARAMS']['assim_freq']
async_freq  = exper['DA_PARAMS']['async_freq']
cook        = exper['DA_PARAMS']['cook']

dt_window   = DT.timedelta(seconds=assim_freq)

if newstart[0] == True:
    start = DT.datetime(newstart[1],newstart[2],newstart[3],newstart[4],newstart[5],newstart[6])
    cook = 0
else:
    cook_window = DT.timedelta(seconds=exper['DA_PARAMS']['cook'])
    start = DT. datetime(exper['YEAR'],
                         exper['MONTH'],
                         exper['DAY'],
                         exper['HOUR'],
                         exper['MINUTE'],
                         exper['SECOND'])

stop = DT.datetime(stop[0],stop[1],stop[2],stop[3],stop[4],stop[5])

print(("\n#==> run_Exper: start time for experiment is %s" % start.strftime("%Y-%m-%d %H:%M:%S")))
print(("\n#==> run_Exper: stop  time for experiment is %s" % stop.strftime("%Y-%m-%d %H:%M:%S")))

time = []


if cook != 0:
    dt = DT.timedelta(seconds=cook)
    time.append([start, 0])
    time.append([start+dt, 2])
else:
     if newstart[0] == True:
         time.append([start, 0])
         time.append([start+dt_window, 2])
     else:
         time.append([start, 2])
         time.append([start+dt_window, 2])
  
while time[-1][0] < stop:
    t = time[-1][0]
    time.append([t+dt_window,2])

# dont want to do correct ens last time before forecast, so always set last value to assim only
time[-1][1] = 1
    
if fcst > 0:
  dt = DT.timedelta(seconds=fcst) 
  if time[-1][0] != stop:  
    time.append([time[-1][0]+dt, 0])
    print("Error: (Stop-Start) time is not modulo assimilation freq time --> fixing this")
  else:
    time.append([stop+dt, 0])

#-------------------------------------------------------------------------------
#
# Print out what you will run....
   
for n, t in enumerate(time[:-1]):
  now_DT   = t[0]
  later_DT = time[n+1][0]
  if t[1] == -1:
    print(("#======>>> Step %d:  Run pre-cook from %s until %s" % (n,now_DT.strftime("%Y_%m-%d %H:%M:%S"),
                                                                 later_DT.strftime("%Y_%m-%d %H:%M:%S"))))
  if t[1] == 0:
    print(("#======>>> Step %2d:  Run forecast from %s until %s" % (n,now_DT.strftime("%Y_%m-%d %H:%M:%S"),
                                                                 later_DT.strftime("%Y_%m-%d %H:%M:%S"))))
  if t[1] == 1:
    print(("#======>>> Step %2d:  Assimilate at %s, then run forecast until %s" % (n,now_DT.strftime("%Y_%m-%d %H:%M:%S"),
                                                                                later_DT.strftime("%Y_%m-%d %H:%M:%S"))))
  if t[1] == 2:
    print(("#======>>> Step %2d:  Assimilate and run additive noise at %s, then run forecast until %s"  \
                                                                             % (n,now_DT.strftime("%Y_%m-%d %H:%M:%S"),
                                                                                later_DT.strftime("%Y_%m-%d %H:%M:%S"))))

#-------------------------------------------------------------------------------
#
# Set up Hx and Assim_Window parameters

if async_freq > 0:
  overshoot = ((int(assim_freq/async_freq))/2)*async_freq    
else:
  overshoot = 0

if overshoot == 0:
  print(("\n#==> run_Exper: Synchronous DA is requested, local window = %d sec\n" % assim_freq))
else:
  print(("\n#==> run_Exper: Asynchronous DA is requested, overshoot = %d sec\n" % overshoot))

#-------------------------------------------------------------------------------
# 
# Main time loop
#

cpu_total = cpu.time()

for n, t in enumerate(time[:-1]):
  
  now_DT   = t[0]
  later_DT = time[n+1][0]
  now      = (now_DT - start).seconds
  action   = t[1]
 
 # if on last step dont do overshoot
 
  if later_DT == time[n+1][0] == time[-1][0]:
      later    = (later_DT - start).seconds
  else:
      later    = (later_DT - start).seconds + overshoot
    
  print("\n#--------------------------------------------------------------------------")
  print(("\n#==> run_Exper: Now   Model TIME is  %s  seconds " % now))
  print(("\n#==> run_Exper: Later Model TIME is  %s  seconds " % later))

  if t[1] == 0:
    print(("\n#==> run_Exper: Step %2d:  Run forecast from %s until %s" \
           % (n,now_DT.strftime("%Y_%m-%d %H:%M:%S"), later_DT.strftime("%Y_%m-%d %H:%M:%S"))))
  if t[1] == 1:
    print(("\n#==> run_Exper: Step %2d:  Assimilate at %s, then run forecast until %s" 
          % (n,now_DT.strftime("%Y_%m-%d %H:%M:%S"), later_DT.strftime("%Y_%m-%d %H:%M:%S"))))
          
  if t[1] == 2:
    print(("\n#==> run_Exper: Step %2d:  Assimilate and run additive noise at %s, then run forecast until %s"  \
          % (n,now_DT.strftime("%Y_%m-%d %H:%M:%S"), later_DT.strftime("%Y_%m-%d %H:%M:%S"))))
#-------------------------------------------------------------------------------
#     
# Assimilate/if requested run additive noise
#
  if action >= 1:  # assimilate observations.....
     
    c0 = cpu.time()
    
    print(("\n#==> run_Exper:  Assimilation is being called, the time:  %s" % now_DT.strftime("%Y-%m-%d %H:%M:%S")))
      
    cmd = "run_filter.py --exper %s --time %s" % (options.exp, now_DT.strftime("%Y,%m,%d,%H,%M,%S"))
    
    if async_freq > 0:  cmd = "%s --freq %d" % (cmd, async_freq)

    if HF[0] == True:  
      cmd = "%s --HF %d" % (cmd, assim_count%HF[1])
      assim_count = assim_count + 1

    cwd = os.getcwd()
    newcmd = os.path.join(cwd, cmd)

    print("\n"+newcmd+"\n")
    
    if RunIt:  
      run_unix_cmd(newcmd)
   #  print(("\nSTDOUT: %s " % master_output)) 
   #  print(("\nSTDERR: %s " % master_error))
    
    cpu_enkf = cpu_enkf + cpu.time() - c0

  if action == 2 and exper['DA_PARAMS']['additive_noise'][0]:  # run additive noise
  
    c0 = cpu.time()
    
    print(("\n#==> run_Exper:  Additive noise is being called, the time: %s" % now_DT.strftime("%Y-%m-%d %H:%M:%S")))
    print(("\n#==> run_Exper:  Additive noise is being called, the time: %s" % now_DT.strftime("%Y-%m-%d %H:%M:%S")))
    
    if exper['DA_PARAMS']['additive_noise'][1] <= 1:
        print("\n#==> run_Exper:  Additive noise is based on composite reflectivity")
        cmd = "ens.py -e %s --crefperts -t %s --write" % (options.exp, now_DT.strftime("%Y,%m,%d,%H,%M,%S"))
    if exper['DA_PARAMS']['additive_noise'][1] == 2:
        print("\n#==> run_Exper:  Additive noise is based on LETKF adaptive inflation field")
        cmd = "ens.py -e %s --ainflperts -t %s --write" % (options.exp, now_DT.strftime("%Y,%m,%d,%H,%M,%S"))

    cwd = os.getcwd()
    newcmd = os.path.join(cwd, cmd)

    print(("\n\n"+newcmd+"\n"))

    if RunIt:  
      os.system(newcmd)
   #  master_output, master_error = run_unix_cmd(newcmd)
   #  print(("\nSTDOUT: %s " % master_output)) 
   #  print(("\nSTDERR: %s " % master_error))

    cpu_correct = cpu_correct + cpu.time() - c0

#-------------------------------------------------------------------------------
#     
# Integrate ensemble members to next observation time.
#
  if action >= 0:  # run model
  
    c0 = cpu.time()

    print(("\n#==> run_Exper: Integrating ensemble members from date and time: %s" % (now_DT.strftime("%Y-%m-%d_%H:%M:%S"))))

    runtime = (later_DT - now_DT).seconds + overshoot
    cmd = "run_fcst.py -e %s --run_time %d -t %s" % (options.exp, runtime, now_DT.strftime("%Y,%m,%d,%H,%M,%S"))
    cwd = os.getcwd()
    newcmd = os.path.join(cwd, cmd)

    print(("\n\n"+newcmd+"\n"))
    
    if RunIt:  
      master_output, master_error = run_unix_cmd(newcmd)
      print(("\nSTDOUT: %s " % master_output)) 
      print(("\nSTDERR: %s " % master_error))

    print(("\n#==> run_Exper: Integrated ensemble members to time: %s\n" % ( later_DT.strftime("%Y-%m-%d_%H:%M:%S"))))
    
    cpu_model = cpu_model + cpu.time() - c0

cpu_total = cpu.time() - cpu_total

#-------------------------------------------------------------------------------

print(("\n#==> run_Exper:  Elapsed WALL CLOCK TIME FOR TOTAL EXPERIMENT: %f \n" % ( cpu_total )))
print(("\n#==> run_Exper:  Elapsed WALL CLOCK TIME FOR FILTER:           %f \n" % ( cpu_enkf )))
print(("\n#==> run_Exper:  Elapsed WALL CLOCK TIME FOR MODEL RUNS:       %f \n" % ( cpu_model )))
print(("\n#==> run_Exper:  Elapsed WALL CLOCK TIME FOR ADDITIVE NOISE:   %f \n" % ( cpu_correct )))

#-------------------------------------------------------------------------------
