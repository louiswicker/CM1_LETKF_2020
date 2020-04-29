#!/usr/bin/env python

import os

cmd = "f2py --fcompiler='gnu95' --f90flags='-O2' -c -m recursive2d recursive2d.f90"
os.system(cmd)
