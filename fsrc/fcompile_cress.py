#!/usr/bin/env python

import os

cmd = "f2py --fcompiler='gnu95' --f90flags='-O3' -c -m cressman cressman.f90"
os.system(cmd)
#cmd = "f2py --fcompiler='gnu95' --f90flags='-O3' -c -m -m uhcalc uhcalc.f90"
#os.system(cmd)
