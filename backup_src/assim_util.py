import sys
import os
import numpy as N
from time import time as timer
import netCDF4 as ncdf

import datetime as DT

# Parameters for the data assimilation

_missing_value = 9900000000.
_time_units    = 'seconds since 1970-01-01 00:00:00'
_calendar      = 'standard'
_stringlen     = 8
_datelen       = 19


# ==================================================================================================================================
class myTimer():

  def __init__(self, start = 0.0, name = 'nothing', minutes = False, args=[], kwargs={}):
    self.name    = name
    self.last    = 0.0
    self.total   = 0.0
    self.minutes = minutes

  def start(self):
    self.last = timer()
      
  def stop(self):
    self.total = self.total + timer() - self.last

  def printit(self,string=None):
    print "\n  ---------------------------------------\n"
    if string == None: string = ""
    total = self.total
    if self.minutes:  
      total = round(total/60., 3) 
      print "  %s ... Time in minutes for %s: %f" % (string, self.name, total)
    else:
      total = round(total, 3) 
      print "  %s ... Time in seconds for %s: %f " % (string, self.name, total)

    print "\n  ---------------------------------------\n"


# ==================================================================================================================================
def wrap_letkf_args(args):
  """ Simple wrapper around last_dimension_stats that works with a single
      argument so that it works with pool.map.  It expects a tuple with
      (start, stop, npy_file) as its input.
  """
  return letkf(*args)
  
# ==================================================================================================================================
def print_state_diagnostics(mdata, store=False, header=False):

  global uavg, vavg, wavg, tavg, davg, qavg    # use global storage so that we can simply store this stuff inside

  if header:
    print "          Z |           UA                    VA                  WA                  TH                  DBZ                QR"
    print "                    mean  std             mean  std           mean  std           mean  std           mean   std         mean  std\n"

  if store:
    uavg = N.zeros((100,2))
    vavg = N.zeros((100,2))
    wavg = N.zeros((100,2))
    tavg = N.zeros((100,2))
    davg = N.zeros((100,2))
    qavg = N.zeros((100,2))

    for zi in N.arange(mdata.nz):
      uavg[zi,0], uavg[zi,1] = N.average(mdata["UA"].data[:-1,zi,:,:]),  N.std(mdata["UA"].data[:-1,zi,:,:])
      vavg[zi,0], vavg[zi,1] = N.average(mdata["VA"].data[:-1,zi,:,:]),  N.std(mdata["VA"].data[:-1,zi,:,:])
      wavg[zi,0], wavg[zi,1] = N.average(mdata["WA"].data[:-1,zi,:,:]),  N.std(mdata["WA"].data[:-1,zi,:,:])
      tavg[zi,0], tavg[zi,1] = N.average(mdata["TH"].data[:-1,zi,:,:]),  N.std(mdata["TH"].data[:-1,zi,:,:])
      davg[zi,0], davg[zi,1] = N.average(mdata["DBZ"].data[:-1,zi,:,:]), N.std(mdata["DBZ"].data[:-1,zi,:,:])
      qavg[zi,0], qavg[zi,1] = N.average(mdata["QR"].data[:-1,zi,:,:]),  N.std(mdata["QR"].data[:-1,zi,:,:])
    return
        
  for zi in N.arange(mdata.nz):
    print " Prior       |  %8.5f  %8.5f  %10.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f" % \
        (uavg[zi,0], uavg[zi,1], vavg[zi,0], vavg[zi,1],  wavg[zi,0], wavg[zi,1], tavg[zi,0], tavg[zi,1], \
         davg[zi,0], davg[zi,1], 1000.*qavg[zi,0], 1000.*qavg[zi,1])

    print " %12.2f|  %8.5f  %8.5f  %10.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f \n" % \
                         (mdata.zc.data[zi],   \
                         N.average(mdata["UA"].data[:-1,zi,:,:]),  N.std(mdata["UA"].data[:-1,zi,:,:]), \
                         N.average(mdata["VA"].data[:-1,zi,:,:]),  N.std(mdata["VA"].data[:-1,zi,:,:]), \
                         N.average(mdata["WA"].data[:-1,zi,:,:]),  N.std(mdata["WA"].data[:-1,zi,:,:]), \
                         N.average(mdata["TH"].data[:-1,zi,:,:]),  N.std(mdata["TH"].data[:-1,zi,:,:]), \
                         N.average(mdata["DBZ"].data[:-1,zi,:,:]),N.std(mdata["DBZ"].data[:-1,zi,:,:]),
                   1000.*N.average(mdata["QR"].data[:-1,zi,:,:]), 1000.*N.std(mdata["QR"].data[:-1,zi,:,:]))

#####################################################################################################
def rename_weight_file(DT):
        
# create the fileput filename and create new netCDF4 file

  filename = "%s_%s%s" % ("Weights", DT.strftime("%Y-%m-%d_%H:%M:%S"), ".nc" )
  print("moving ./%s to ./%s" % ("weights.nc", filename))
  os.system("mv ./%s ./%s" % ("weights.nc", filename))

  return

#####################################################################################################
def write_inflation_file(array, xc, yc, zc, DT):
        
# create the fileput filename and create new netCDF4 file

  filename = "%s_%s%s" % ("Inflation", DT.strftime("%Y-%m-%d_%H:%M:%S"), ".nc" )

  print " Writing %s as the adaptive inflation file..." % (filename)
    
  rootgroup = ncdf.Dataset(filename, 'w', format='NETCDF4')
      
# Create dimensions

  shape = array.shape
  
# rootgroup.createDimension('nvar', shape[0])
  rootgroup.createDimension('nz',   shape[0])
  rootgroup.createDimension('ny',   shape[1])
  rootgroup.createDimension('nx',   shape[2])
  rootgroup.createDimension('stringlen', _stringlen)
  rootgroup.createDimension('datelen', _datelen)
  
# Write some attributes

  rootgroup.time_units = _time_units
  rootgroup.calendar   = _calendar
  rootgroup.stringlen  = "%d" % (_stringlen)
  rootgroup.datelen    = "%d" % (_datelen)

# Create variables

# V_type  = rootgroup.createVariable('inflation', 'f8', ('nvar', 'nz', 'ny', 'nx'), zlib=True, shuffle=True )    
  V_type  = rootgroup.createVariable('inflation', 'f8', ('nz', 'ny', 'nx'), zlib=True, shuffle=True )    
  V_dates = rootgroup.createVariable('date', 'S1', ('datelen'), zlib=True, shuffle=True)
# V_nv    = rootgroup.createVariable('VAR', 'i4', ('nvar'), zlib=True, shuffle=True)
  V_xc    = rootgroup.createVariable('XC', 'f4', ('nx'), zlib=True, shuffle=True)
  V_yc    = rootgroup.createVariable('YC', 'f4', ('ny'), zlib=True, shuffle=True)
  V_zc    = rootgroup.createVariable('ZC', 'f4', ('nz'), zlib=True, shuffle=True)

# Write variables

  rootgroup.variables['date'][:] = ncdf.stringtoarr(DT.strftime("%Y-%m-%d_%H:%M:%S"), _datelen)
  
  rootgroup.variables['inflation'][:] = array[:]
  rootgroup.variables['XC'][:]        = xc[:]
  rootgroup.variables['YC'][:]        = yc[:]
  rootgroup.variables['ZC'][:]        = zc[:]
# rootgroup.variables['VAR'][:]       = N.arange(shape[0])
  
  rootgroup.sync()
  rootgroup.close()
  
  return

#####################################################################################################
def read_inflation_file(DT=None, file=None):
        
# create the fileput filename and create new netCDF4 file

  if file != None:
    filename = file
  elif DT != None:
    filename = "%s_%s%s" % ("Inflation", DT.strftime("%Y-%m-%d_%H:%M:%S"), ".nc" )
  else:
    print "Read_inflation_file:  MUST SPECIFY FILENAME OR DateTIME object"
    print "====>Exiting!!!"
    return -1
    
  rootgroup = ncdf.Dataset(filename, 'r', format='NETCDF4')
  
  return N.asfortranarray(rootgroup.variables['inflation'][:])

# ==================================================================================================================================
def write_ob_table():
#
## Define observation types and how they are to be used during the assimilation.

  ObTableFile = "ob_table.txt"
  
  ObTable = """17    ! number of observation types listed below
  'U10m'     'm/s'    1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'V10m'     'm/s'    1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'TEMP2m'   'K'      1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'DEWPT2m'  'K'      1    2.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'QV2m'     'g/g'    1    0.001   0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'U'        'm/s'    1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'V'        'm/s'    1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'W'        'm/s'    1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'THETA'    'K'      1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'P'        'mb'     1    1.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'QV'       'g/g'    1    0.001   0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'QR'       'g/g'    1    0.001   0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'ZDR'      'dB'     1    0.6     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'KDP'      'deg/km' 1    4.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'VR'       'm/s'    1    2.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'REFL'     'dBZ'    1    5.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  'FLSH2D'   'fl/s'   1    2.0     0.0    1  1  1  0  0  1  1  1  1  1  1  1  1   1   1   1   1  1   1   0   0   0    0   0    0    0    0   0    1   1   1   1   1   1   1   1   1   1   1   1   1   1    1
  ---------  -----    ---  ------  ----   -- -- -- -- -- -- -- -- -- -- -- -- --- --- --- --- -- --- --- --- --- ---- --- ---- ---- ---- --- ---- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  ---
  ob. type   units    use  random  bias   U  V  W  PI KM TH QV QC QR QI QS QH QIR QGL QGM QGH QF QIP QHL RTC RTI RTIR RTS RTGL RTGM RTGH RTF RTIP CCW CRW CIP CSW CHW CHL VHW VHL QHW QSW ZRW ZSW ZHW ZHL  QHLW
                          error   error  update model field listed above?  (1=yes, 0=no)
                          """
  ofile = open(ObTableFile, 'w+')
  ofile.write(ObTable)
  ofile.close()
  
  
#####################################################################################################
# Main program
if __name__ == "__main__":

  
  print "Testing I/O of inflation array"
  
  dt1= DT.datetime(2003,5,8,22,1,2)

  inflate1 = N.array(N.random.random((50,100,120)), order='F', dtype=N.float64)

  x, y, z = N.arange(inflate1.shape[2]), N.arange(inflate1.shape[1]), N.arange(inflate1.shape[0])

  write_inflation_file(inflate1,x,y,z,dt1)

  inflate2 = read_inflation_file(DT=dt1)

  if N.any(inflate1-inflate2 != 0.0):
    print "====> problem!, inflation arrays do not match!!"
  else:
    print "====> good to go as the inflation arrays match!!"

