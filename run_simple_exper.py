#!/usr/bin/env python

import os

print("\nStarting run_simple_exper.py")

#--------------------------------------------------------------------------
# Initialize the run, and add perturbations to the backround and 3D fields

os.system("create_run_letkf.py ")
os.system("run_fcst.py -e RUN_LETKF/RUN_LETKF.exp -i ")
os.system("ens.py -e RUN_LETKF/RUN_LETKF.exp --init0 -t 2003,5,8,20,40,0 --write ")

#--------------------------------------------------------------------------
# this is the cook period....

os.system("/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 1200 -t 2003,05,08,20,40,00")

#--------------------------------------------------------------------------
# Now loop through the cycling at 5 min intervals

times = ['2003,05,08,21,00,00', '2003,05,08,21,05,00', '2003,05,08,21,10,00', '2003,05,08,21,15,00', \
         '2003,05,08,21,20,00', '2003,05,08,21,25,00', '2003,05,08,21,30,00', '2003,05,08,21,35,00', \
         '2003,05,08,21,40,00', '2003,05,08,21,45,00', '2003,05,08,21,50,00', '2003,05,08,21,55,00', '2003,05,08,22,00,00']

for time in times:

    os.system("/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp -t %s --freq -300" % (time))

    if time != times[-1]:  # dont the crefs or last 5 min forecast - do the pure forecast at the end

#       os.system("/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t %s --write" % (time))

        os.system("/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t %s" % (time))

#--------------------------------------------------------------------------
# Make a 30 minute forcast

    os.system("/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 1800 -t 2003,5,8,22,00,00")

#--------------------------------------------------------------------------
# Completed (hopefully) simple synchronous experiment

#--------------------------------------------------------------------------
# Make a few plots every 10 minutes

for time in times[::2]:
    os.system("/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp -t %s -v W --plot9" % (time))
    os.system("/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp -t %s -v WZ --plot9" % (time))
    os.system("/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp -t %s -v DBZ --plot8" % (time))

#--------------------------------------------------------------------------
# Creat DA diagnostics 

os.system("python DBZ_CR.py -d RUN_LETKF -t DBZ_CR --noshow")
os.system("python DBZ_INV.py -d RUN_LETKF -t DBZ_INV --noshow")
os.system("python VR_CR.py -d RUN_LETKF -t VR_CR --noshow")
os.system("python VR_INV.py -d RUN_LETKF -t VR_INV --noshow")

print("\nEnded run_simple_exper.py")

