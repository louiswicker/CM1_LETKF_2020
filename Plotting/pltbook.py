import numpy as N
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

def nice_mxmnintvl( dmin, dmax, outside=False, max_steps=20, cint=None, sym=False):

    """ Description: Given min and max values of a data domain and the maximum
                     number of steps desired, determines "nice" values of
                     for endpoints and spacing to create a series of steps
                     through the data domainp. A flag controls whether the max
                     and min are inside or outside the data range.

        In Args: float   dmin             the minimum value of the domain
                 float   dmax       the maximum value of the domain
                 int     max_steps  the maximum number of steps desired
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
        return None, None, N.zero(1)

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
    #print am1
    #print ax1
    #print cints

    if cint == None or cint == 0.0:
        try:
            index = N.where(cints < max_steps)[0][0]
            return am1[index], ax1[index], N.linspace(am1[index], ax1[index], cints[index])
        except IndexError:
            return None, None, N.zero(1)
    else:
        return am1, ax1, cint
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
