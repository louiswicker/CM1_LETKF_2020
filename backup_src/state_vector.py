#-----------------------------------------------------------------------------------------
# List of supported state vectors for full state data assimilation reads
#      ENS.PY interates over the states if no state vector is supplied to try and match
#      microphysics attribute in netCDF file --> therefore the names of the state declarations
#      need to match the MICROPHYS global attribute found the netCDF COMMAS files
#
#      At the bottom, there is a special HxF state vector for only reading in variables
#      and coordinate information needed to create observation priors/posteriors.  This
#      cuts down on the I/O time needed for asynchronous DA

#-----------------------------------------------------------------------------------------
# BEGIN STATE DECLARATION

morrison = { 
         "nxyz3d": 12,
         "nxy2d":  0,
         "nxz2d":  0,
         "nyz2d":  0,
         "xyz3d":  ["U", "V", "W", "DBZ", "PI0", "TH", "QV", "QC", "QR", "QI", "QS", "QG"],
         "coords": ["xc", "yc", "zc", "xe", "ye", "ze"],
         "U"     : {"name": "ua",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "V"     : {"name": "va",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "W"     : {"name": "wa",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "DBZ"   : {"name": "dbz",    "ndims": 3, "posdef": 1, "inflation": 0, "writeback": True},
         "PI0"   : {"name": "pi0",    "ndims": 3, "posdef": 0, "inflation": 0, "writeback": False},
         "TH"    : {"name": "theta",  "ndims": 3, "posdef": 0, "inflation": 2, "writeback": True},
         "QV"    : {"name": "qv",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QC"    : {"name": "qc",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QR"    : {"name": "qr",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QI"    : {"name": "qi",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QS"    : {"name": "qs",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QG"    : {"name": "qg",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "xc"    : {"name": "xh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "yc"    : {"name": "yh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "zc"    : {"name": "zh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "xe"    : {"name": "xf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "ye"    : {"name": "yf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "ze"    : {"name": "zf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
        }

# END STATE DECLARATION
#-----------------------------------------------------------------------------------------
# BEGIN STATE DECLARATION

zvdLFO = { 
         "nxyz3d": 12,
         "nxy2d":  0,
         "nxz2d":  0,
         "nyz2d":  0,
         "xyz3d":  ["U", "V", "W", "DBZ", "PI0", "TH", "QV", "QC", "QR", "QI", "QS", "QG"],
         "coords": ["xc", "yc", "zc", "xe", "ye", "ze"],
         "U"     : {"name": "ua",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "V"     : {"name": "va",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "W"     : {"name": "wa",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "DBZ"   : {"name": "dbz",    "ndims": 3, "posdef": 1, "inflation": 0, "writeback": True},
         "PI0"   : {"name": "pi0",    "ndims": 3, "posdef": 0, "inflation": 0, "writeback": False},
         "TH"    : {"name": "theta",  "ndims": 3, "posdef": 0, "inflation": 2, "writeback": True},
         "QV"    : {"name": "qv",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QC"    : {"name": "qc",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QR"    : {"name": "qr",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QI"    : {"name": "qi",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QS"    : {"name": "qs",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QG"    : {"name": "qg",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "xc"    : {"name": "xh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "yc"    : {"name": "yh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "zc"    : {"name": "zh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "xe"    : {"name": "xf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "ye"    : {"name": "yf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "ze"    : {"name": "zf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
        }

# END STATE DECLARATION
#-----------------------------------------------------------------------------------------
# BEGIN STATE DECLARATION

zvdh = { 
         "nxyz3d": 13,
         "nxy2d":  0,
         "nxz2d":  0,
         "nyz2d":  0,
         "xyz3d":  ["U", "V", "W", "DBZ", "PI0", "TH", "QV", "QC", "QR", "QI", "QS", "QG", "QH"],
         "coords": ["xc", "yc", "zc", "xe", "ye", "ze"],
         "U"     : {"name": "ua",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "V"     : {"name": "va",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "W"     : {"name": "wa",     "ndims": 3, "posdef": 0, "inflation": 0, "writeback": True},
         "DBZ"   : {"name": "dbz",    "ndims": 3, "posdef": 1, "inflation": 0, "writeback": True},
         "PI0"   : {"name": "pi0",    "ndims": 3, "posdef": 0, "inflation": 0, "writeback": False},
         "TH"    : {"name": "theta",  "ndims": 3, "posdef": 0, "inflation": 2, "writeback": True},
         "QV"    : {"name": "qv",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QC"    : {"name": "qc",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QR"    : {"name": "qr",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QI"    : {"name": "qi",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QS"    : {"name": "qs",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QG"    : {"name": "qg",     "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "QH"    : {"name": "qhl",    "ndims": 3, "posdef": 1, "inflation": 3, "writeback": True},
         "xc"    : {"name": "xh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "yc"    : {"name": "yh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "zc"    : {"name": "zh",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "xe"    : {"name": "xf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "ye"    : {"name": "yf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
         "ze"    : {"name": "zf",     "ndims": 1, "posdef": 0, "inflation": 0, "writeback": True},
        }

# END STATE DECLARATION
#-----------------------------------------------------------------------------------------
# BEGIN STATE DECLARATION for reading in ensemble data to compute radar forward operators

Hxf = { 
							"nxyz3d": 5,
							"nxy2d":  0,
							"nxz2d":  0,
							"nyz2d":  0,
							"xyz3d":  ["U", "V", "W", "DBZ", "DEN"],
							"coords": ["xc", "yc", "zc", "xe", "ye", "ze"],
							"U"     : {"name": "ua",     "ndims": 3, "posdef": 0, "inflation": 0},
							"V"     : {"name": "va",     "ndims": 3, "posdef": 0, "inflation": 0},
							"W"     : {"name": "wa",     "ndims": 3, "posdef": 0, "inflation": 0},
							"DBZ"   : {"name": "dbz",    "ndims": 3, "posdef": 1, "inflation": 0},
							"DEN"   : {"name": "rho0",   "ndims": 3, "posdef": 1, "inflation": 0},
							"xc"    : {"name": "xh",     "ndims": 1, "posdef": 0, "inflation": 0},
							"yc"    : {"name": "yh",     "ndims": 1, "posdef": 0, "inflation": 0},
							"zc"    : {"name": "zh",     "ndims": 1, "posdef": 0, "inflation": 0},
							"xe"    : {"name": "xf",     "ndims": 1, "posdef": 0, "inflation": 0},
							"ye"    : {"name": "yf",     "ndims": 1, "posdef": 0, "inflation": 0},
							"ze"    : {"name": "zf",     "ndims": 1, "posdef": 0, "inflation": 0},
        }

# END OF HXF STATE DECLARATION

#-----------------------------------------------------------------------------------------
# BEGIN STATE DECLARATION for reading in ensemble data to check initialization

Init = { 
							"nxyz3d": 8,
							"nxy2d":  0,
							"nxz2d":  0,
							"nyz2d":  0,
							"xyz3d":  ["U", "V", "W", "DBZ", "DEN", "PI0", "TH0", "QV0"],
							"coords": ["xc", "yc", "zc", "xe", "ye", "ze"],
							"U"     : {"name": "ua",     "ndims": 3, "posdef": 0, "inflation": 0},
							"V"     : {"name": "va",     "ndims": 3, "posdef": 0, "inflation": 0},
							"W"     : {"name": "wa",     "ndims": 3, "posdef": 0, "inflation": 0},
							"DBZ"   : {"name": "dbz",    "ndims": 3, "posdef": 1, "inflation": 0},
							"DEN"   : {"name": "rho0",   "ndims": 3, "posdef": 1, "inflation": 0},
							"PI0"   : {"name": "pi0",    "ndims": 3, "posdef": 0, "inflation": 0},
							"TH0"   : {"name": "th0",    "ndims": 3, "posdef": 0, "inflation": 0},
							"QV0"   : {"name": "qv0",    "ndims": 3, "posdef": 1, "inflation": 0},
							"xc"    : {"name": "xh",     "ndims": 1, "posdef": 0, "inflation": 0},
							"yc"    : {"name": "yh",     "ndims": 1, "posdef": 0, "inflation": 0},
							"zc"    : {"name": "zh",     "ndims": 1, "posdef": 0, "inflation": 0},
							"xe"    : {"name": "xf",     "ndims": 1, "posdef": 0, "inflation": 0},
							"ye"    : {"name": "yf",     "ndims": 1, "posdef": 0, "inflation": 0},
							"ze"    : {"name": "zf",     "ndims": 1, "posdef": 0, "inflation": 0},
        }

# END OF HXF STATE DECLARATION
