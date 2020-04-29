from __future__ import division
import numpy as N
import datetime
import netCDF4

hscale = 1.0 
debug  = False

#===============================================================================
def GetModelTime(file_obj, DateTime, index=False, closest=False):
    """ file_obj is a netCDF data set pointer, and DateTime is a tuple of (2003,5,8,22,10,0)

        GetModelTime can either return the corresponding model time in seconds from a DT tuple,
        or can return the index into the time coordinate by setting the index key to True """
   
    DT = datetime.datetime(int(DateTime[0]),int(DateTime[1]),int(DateTime[2]),int(DateTime[3]),int(DateTime[4]),int(DateTime[5]))

    tarray    = file_obj.variables['TIME'][:]
    init_time = file_obj.COARDS.split('since')[1][1:]                                      # strip out COARDS info
    init      = datetime.datetime.strptime(init_time, '%Y-%m-%d %H:%M:%S')                 # create a datetime object
    DateTimeArray  = [init + datetime.timedelta(seconds=int(s)) for s in tarray ]  # this is a list!

# Now find the correct time_index to read from using the datetime objects..

    try:
        coards = file_obj.COARDS.split('since')[1][1:]
        coards = file_obj.COARDS[:]
    except:
        return 'Error Parsing COARDS'

    if debug:
        print "\n COARDS:  ",coards
        print "Input date and time:  ",DT,"\n"
        print "Model date and time:  ",DateTimeArray,"\n"

    time = int(round(netCDF4.date2num(DT, coards)))

    if index:
        return GetTimeIndex(fig_obj, time, closest=closest)
    else:
        return time

    return None

#===============================================================================
def GetTimeIndex(file_obj, time, closest=False):
    """ Returns the index of the time coordinate for a given model time"""

    tarray    = file_obj.variables['TIME'][:]

    index = N.where(tarray >= time)
    if closest:
            index2 = N.where(tarray < time)
            if N.size(index2) !=0 and N.size(index) !=0:
                diff1 = tarray[index[0][0]] - time
                diff2 = time - tarray[index2[0][-1]]
                if diff2 < diff1:
                    return index2[0][-1]
    else:
       if N.size(index) !=0:
           return index[0][0]

#===============================================================================
def GetFileDatetime(file_obj, time, string=False):
    """ Returns a string of the date time from model time index"""

    try:
        init_time = file_obj.COARDS.split('since')[1][1:]
    except:
        return 'Error Parsing COARDS'

    init = datetime.datetime.strptime(init_time, '%Y-%m-%d %H:%M:%S')

    timeindex = GetTimeIndex(file_obj, time)

    try:
        offset = file_obj.variables['TIME'][timeindex]
    except:
        return 'Error Parsing netcdf file time variable'

    valid_time = init + datetime.timedelta(seconds=int(offset))

    if string:
        return valid_time.strftime('%Y-%m-%d %H:%M:%S')
    else:
        return valid_time
 
#===============================================================================
def ComputeWZ(x, y, u, v):
    """ Returns an array of vertical vorticity on A-grid 
        U & V can be 3D arrays, and a 3D vertical vorticity
        volume will be returned at the cell centers """

    dx   = x[1]-x[0]
    dy   = y[1]-y[0]
    nx   = x.shape[0]
    ny   = y.shape[0]

    if len(u.shape) == 3:  # 3D volume of vorticity
        wz = N.zeros((u.shape[0],nx+1,ny+1))
        wz[:,1:nx,1:ny] = (v[:,1:nx,1:ny+1] - v[:,0:nx-1,1:ny+1]) / dx \
                        - (u[:,1:nx,1:ny]   - u[:,1:nx,0:ny-1])   / dy
        return 0.25*(wz[:,0:nx,0:ny] + wz[:,1:nx+1,0:ny] + wz[:,0:nx,1:ny+1] + wz[:,1:nx+1,1:ny+1])

    else:
        wz = N.zeros((nx+1,ny+1))
        wz[1:nx,1:ny] = (v[1:nx,1:ny+1] - v[0:nx-1,1:ny+1]) / dx \
                      - (u[1:nx,1:ny]   - u[1:nx,  0:ny-1]) / dy
        return 0.25*(wz[0:nx,0:ny] + wz[1:nx+1,0:ny] + wz[0:nx,1:ny+1] + wz[1:nx+1,1:ny+1])


# if we get this far, its an Error...
    return N.nan
    
#===============================================================================
def GetVR(x, y, z, rlat, rlon, glat, glon, u, v, w, dbz = None, min_dis=500.):
    """ Returns a 2-D array of radial velocity """

    rdx, rdy = dll_2_dxy(glat, rlat, glon, rlon, degrees=True)

    dx   = x[1]-x[0]
    dy   = y[1]-y[0]
    grid = N.mgrid[x[0]:x[-1]+dx:dx, y[0]:y[-1]+dx:dy]
    x2d  = grid[0] - rdx 
    y2d  = grid[1] - rdy 

    r = N.sqrt(x2d**2 + y2d**2 + z**2)

    return N.where(r >= min_dis, (u*x2d + v*y2d + w*z)/r, 0.0)
    
#===============================================================================
def GetWindSpeed(u,v, w = None):
    """ Returns the array of wind speeds. """
    
    if w == None:
        return (u**2 + v**2)**0.5
    else:
        return (u**2 + v**2 + w**2)**0.5
    
#===============================================================================
def UVW2aGrid(U_edge_points, V_edge_points, W_edge_points = None):
    """ Function that takes both the U edge points and the V edge points
        and returns two arrays that have averaged the edge points and are
        centered on the center point of each grid box. """
    
    u_tmp = (U_edge_points[:,:-1] + U_edge_points[:,1:]) / 2.0
    v_tmp = (V_edge_points[:-1,:] + V_edge_points[1:,:]) / 2.0
 
    if W_edge_points != None:    
      w_tmp = (W_edge_points[0,:,:] + W_edge_points[1,:,:]) / 2.0

    if W_edge_points == None:    
      return (u_tmp, v_tmp)
    else:
      return (u_tmp, v_tmp, w_tmp)
    
#===============================================================================
def GetEnsembleFilenames(file_name, members, i=0, j=0):
    """ Function that takes a single filename and returns the correct
        filename of the ensemble. """
        
    tmp_name = file_name.split('.')
    file_name = ''
    for piece in tmp_name[:-2]:
        file_name += piece
        file_name += '.'
    file_name += '%03i.%s' % (members[i][j], tmp_name[-1])
    print file_name
    
    return file_name

#===============================================================================
def GetDensity(file_obj):
    """     """
    piinit = file_obj.variables['PIINIT'][:]
    thinit = file_obj.variables['THINIT'][:]

    return piinit**2.5 * 1.0e5 / (287.04 * thinit)

#===============================================================================
def GetCoords(file_obj, XC='XC', YC='YC', XE='XE', YE='YE', ZC='ZC', ZE='ZE', tindex=None):
    """     """
    
    xc = file_obj.variables[XC][:]
    yc = file_obj.variables[YC][:]
    zc = file_obj.variables[ZC][:]
    xe = file_obj.variables[XE][:]
    ye = file_obj.variables[YE][:]
    ze = file_obj.variables[ZE][:]

    if tindex != None:  # get offset for moving grid
        x_sw  = file_obj.variables['XG_POS'][tindex]
        y_sw  = file_obj.variables['YG_POS'][tindex]
        xc = xc + x_sw
        xe = xe + x_sw
        yc = yc + y_sw
        ye = ye + y_sw

    return ((xc,yc,zc),(xe,ye,ze))

#===============================================================================
def calcUH(w, vort, zc, ze, zbot, ztop):

    count = 0
    nx=w.shape[2]
    ny=w.shape[1]

    k1 = N.where(zc > zbot)[0][0] 
    k2 = N.where(zc > ztop)[0][0]
    
    print "CalcUH:  nx,ny:  ", nx,ny
    print "CalcUH:  Layer is:  ", ze[k1], ze[k2]
    print "CalcUH:  Max/Min W:     ", w.max(), w.min()
    print "CalcUH:  Max/Min VORT:  ", vort.max(), vort.min()

    UH = -1.0 * N.ones((ny, nx))

    for j in range(1,ny-1):
        for i in range(1,nx-1):
            for k in range(k1, k2+1):
                if (vort[k,j,i] > 0.0 and (w[k,j,i]+w[k+1,j,i]) > 0.0):
                    UH[j,i] = UH[j,i] + 0.5*(w[k,j,i]+w[k+1,j,i])*vort[k,j,i]*(ze[k+1]-zc[k])

#           index_w = (w[k,:,:] > 0.0)
#           index_z = (vort[k,:,:] > 0.0) 
#           index   = index_w & index_z
#           print k, index_w.shape, index_z.shape, index.shape
#           print k, w[k,index].shape
#           UH[:,:] = UH[:,:] + .5*(w[k,index]+w[k+1,index]).reshape(ny,nx)*vort[k,index].reshape(ny,nx)

    print "CalcUH:  Max value:  ", UH.max()
    return UH

#===============================================================================
def get_loc(x0, xc, radius):
  """get_loc returns a range of indices for a coordinate system, a pt, and a radius
  """

  indices = N.where(x-radius <= xc <= x+radius)

  if N.size(indices[0]) == 0:
    return -1, 0
  else:
    i0 = indices[0][0]
    i1 = indices[0][-1]
    return i0, i1

#===============================================================================
def interp_weights(x, xc, extrapolate=False):
  """interp_weights returs the linear interpolation weights for a given
     ascending array of coordinates.

     x = location to be interpolated to
     xc = locations of grid
     extrapolate:  what to do at the boundary?
                   for now, if at edge, set return values to missing

     OUTPUTS:  i0, i1, dx0, dx1 locations and weights for the interpolation
  """

  indices = N.where(xc <= x)

  if N.size(indices[0]) == 0:
    return -1, -1, None, None, None
  else:
    i0 = indices[0][-1]
    if i0 == N.size(xc)-1:
      return -1, -1, None, None, None
    else:
      dx  = xc[i0+1] - xc[i0]
      dx0 = xc[i0+1] - x
      dx1 = dx-dx0
      return i0, i0+1, dx0, dx1, dx

#===============================================================================
def nearlyequal(a, b, sig_digit=None):
    """ Measures the equality (for two floats), in unit of decimal significant 
        figures.  If no sigificant digit is specified, default is 7 digits. """

    if sig_digit == None or sig_digit > 7:
        sig_digit = 7
    if a == b:
        return True
    difference = abs(a - b)
    avg = (a + b)/2
    
    return N.log10(avg / difference) >= sig_digit
    
    
def nice_mxmnintvl( dmin, dmax, outside=True, max_steps=15, cint=None, sym=False):
    """ Description: Given min and max values of a data domain and the maximum
                     number of steps desired, determines "nice" values of 
                     for endpoints and spacing to create a series of steps 
                     through the data domainp. A flag controls whether the max 
                     and min are inside or outside the data range.
  
        In Args: float   dmin 		the minimum value of the domain
                 float   dmax       the maximum value of the domain
                 int     max_steps	the maximum number of steps desired
                 logical outside    controls whether return min/max fall just
                                    outside or just inside the data domainp.
                     if outside: 
                         min_out <= min < min_out + step_size
                                         max_out >= max > max_out - step_size
                     if inside:
                         min_out >= min > min_out - step_size
                                         max_out <= max < max_out + step_size
      
                 float    cint      if specified, the contour interval is set 
                                    to this, and the max/min bounds, based on 
                                    "outside" are returned.

                 logical  sym       if True, set the max/min bounds to be anti-symmetric.
      
      
        Out Args: min_out     a "nice" minimum value
                  max_out     a "nice" maximum value  
                  step_size   a step value such that 
                                     (where n is an integer < max_steps):
                                      min_out + n * step_size == max_out 
                                      with no remainder 
      
        If max==min, or a contour interval cannot be computed, returns "None"
     
        Algorithm mimics the NCAR NCL lib "nice_mxmnintvl"; code adapted from 
        "nicevals.c" however, added the optional "cint" arg to facilitate user 
        specified specific interval.
     
        Lou Wicker, August 2009 """

    table = N.array([1.0,2.0,2.5,4.0,5.0,10.0,20.0,25.0,40.0,50.0,100.0,200.0,
                      250.0,400.0,500.0])

    if nearlyequal(dmax,dmin):
        return None
    
    # Help people like me who can never remember - flip max/min if inputted reversed
    if dmax < dmin:
        amax = dmin
        amin = dmax
    else:
        amax = dmax
        amin = dmin

    if sym:
        smax = max(amax.max(), amin.min())
        amax = smax
        amin = -smax

    d = 10.0**(N.floor(N.log10(amax - amin)) - 2.0)
    if cint == None or cint == 0.0:
        t = table * d
    else:
        t = cint
    if outside:
        am1 = N.floor(amin/t) * t
        ax1 = N.ceil(amax/t)  * t
        cints = (ax1 - am1) / t 
    else:
        am1 = N.ceil(amin/t) * t
        ax1 = N.floor(amax/t)  * t
        cints = (ax1 - am1) / t
    
    # DEBUG LINE BELOW
    # print t, am1, ax1, cints
    
    if cint == None or cint == 0.0:   
        try:
            index = N.where(cints < max_steps)[0][0]
            return am1[index], ax1[index], cints[index]
        except IndexError:
            return None
    else:
        return am1, ax1, cint


def nice_clevels( *args, **kargs):
    """ Extra function to generate the array of contour levels for plotting 
        using "nice_mxmnintvl" code.  Removes an extra step.  Returns 4 args,
        with the 4th the array of contour levels.  The first three values
        are the same as "nice_mxmnintvl". """
    
    try:
        amin, amax, cint = nice_mxmnintvl(*args, **kargs)
        return amin, amax, cint, N.arange(amin, amax+cint, cint) 
    except:
        return None

#===================================================================================================
def compute_az_el(x, y, z, degrees=True):

    rearth  = 6667.e3
    eer     = rearth * 4. / 3.
    rxy     = N.sqrt( x**2 + y**2)

    az = N.nan
    el = N.nan

    if x == 0.0 and y == 0.0:
        az = 0.0

    elif y == 0.0:
        if x > 0.0:
            az = 0.5*N.pi
        else:
            az = 1.5*N.pi

    elif x >= 0.0 and y > 0.0:
        az = N.math.atan(x/y)

    elif x >= 0.0 and y < 0.0:
        az = -N.math.atan(x/abs(y)) + N.pi

    elif x < 0.0 and y < 0.0:
        az = N.math.atan(x/y) + N.pi

    else:
      az = -N.math.atan(abs(x)/y) + 2.0*N.pi

    if rxy == 0.0:
        if z > 0.0:
            el = 0.5*N.pi
        elif z <= 0.0:
            el = -0.5*N.pi
    else:
       el = N.math.atan(z/rxy)

    if degrees:
        return  N.degrees(az), N.degrees(el)
    else:
        return az, el        

#===============================================================================
def dll_2_dxy(lat1, lat2, lon1, lon2, degrees=False, azimuth=False, proj = 'latlon'):

  """dll_2_dxy returns the approximate distance in meters between two lat/lon pairs

     Valid projections: Lambert Conformal
                        Lat - Lon

     INPUTS: Two (lat,lon) pairs in radians, or if degrees==True, degrees (default)

     if lon2 > lon1:  x > 0

     if lat2 > lat1:  y > 0

     OUTPUTS:  DX, DY in meters

     Azimuth formula from http://www.movable-type.co.uk/scripts/latlong.html
  """

  if degrees:
    rlon1 = N.deg2rad(lon1)
    rlon2 = N.deg2rad(lon2)
    rlat1 = N.deg2rad(lat1)
    rlat2 = N.deg2rad(lat2)
  else:
    rlon1 = lon1
    rlon2 = lon2
    rlat1 = lat1
    rlat2 = lat2

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth  = 1000.0 * 6367.0
    x       = rearth * N.cos(0.5*(rlat1+rlat2)) * (rlon2-rlon1)
    y       = rearth * (rlat2-rlat1)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    x, y = p1(lon2, lat2, errchk = True)

  if azimuth:
    ay = N.sin(rlon2-rlon1)*N.cos(rlat2)
    ax = N.cos(rlat1)*N.sin(rlat2)-N.sin(rlat1)*N.cos(rlat2)*N.cos(rlon2-rlon1)
    az = N.degrees(N.arctan2(ay,ax))
    return x, y, az

  return x, y

#===============================================================================
def dxy_2_dll(x, y, lat1, lon1, degrees=True, proj = 'latlon'):

  """dxy_2_dll returns the approximate lat/lon between an x,y coordinate and
     a reference lat/lon point.

     Valid projections: Lambert Conformal
                        Lat - Lon

     INPUTS:  x,y in meters, lat1, lon1 in radians, or if degrees == True,
              then degrees (default value)

     if x > 0, lon > lon1

     if y > 0, lat > lat1

     OUTPUTS:  lat, lon in radians, or if degrees == True, degrees
               i.e., the input and output units for lat/lon are held the same.
  """

  if degrees:
    rlon1 = N.deg2rad(lon1)
    rlat1 = N.deg2rad(lat1)
  else:
    rlon1 = lon1
    rlat1 = lat1

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth = 1000.0 * 6367.0
    rlat2  = rlat1 + y / rearth
    lon    = N.rad2deg(rlon1 + x / ( rearth * N.cos(0.5*(rlat1+rlat1)) ) )
    lat    = N.rad2deg(rlat2)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    lon, lat = p1(x, y, inverse = True)

  if degrees == False:
    return N.deg2rad(lat), N.deg2rad(lon)
  else:
    return lat, lon

#===================================================================================================        
if __name__ == "__main__":
    print
    print "Testing nice_mxmnintvl code....."
    print
    print "Answers should look like this....."
    print "(-20.0, 70.0, 9.0)"
    print "-----OUTPUT---------"
    print nice_mxmnintvl(-23.1, 70.7, outside=False, cint = 0.0)
    print nice_mxmnintvl(-23.1, 70.7, outside=False)
    print
    print "Testing nice_clevels code....."
    print
    print "Answer should look like this....."
    print "(-20.0, 70.0, 5.0, array([-20., -15., -10.,  -5.,   0.,   5.,  10.,  15.,  20.,  25.,  30."
    print "        35.,  40.,  45.,  50.,  55.,  60.,  65.,  70.]))"
    print
    print "-Output-----------------------------------------------------------------------------------"
    print nice_clevels(-23.1, 70.7, outside=False, cint=5.0)
    print
    print "Answer should look like this....."
    print "(-25.0, 75.0, 5.0, array([-25., -20., -15., -10.,  -5.,   0.,   5.,  10.,  15.,  20.,  25."
    print "        30.,  35.,  40.,  45.,  50.,  55.,  60.,  65.,  70.,  75.]))"
    print
    print "-Output-----------------------------------------------------------------------------------"
    print nice_clevels(-23.1, 70.7, outside=True, cint = 5.0)
