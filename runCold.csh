#!/bin/csh
#
set run_anal = 0
set run_1min = 1
set run_5min = 1

set nens     = 54

set prefix = "PAR_COLD"
set fprefix = "par_2km"
set nthreads = 27

#-----------------------------------------------------------------------------------------
#  Creating the 1 min synchronous analyses

if ( $run_anal == 1 ) then
   echo "creating the 1 min analyses" 
   python pyinit2km_cold.py

   python cold_anal.py -f $fprefix.exp >& $fprefix.out
   mkdir "PAR_1min_Analysis"
   mv $fprefix.* "PAR_1min_Analysis"/.
   mv Prior* "PAR_1min_Analysis"/.
   mv Inflation* "PAR_1min_Analysis"/.
   echo "Ended 1 min cycle"
   date

endif

#-----------------------------------------------------------------------------------------
#  PAR experiment

if ( $run_1min == 1 ) then

   echo "STARTING 1 minute experiment.."
#  python pyinit2km_cold.py

#-----------------------------------------------------------------------------------------

#  echo "Running forecast from 20:30 UTC"

#  python cold_1min.py -f $fprefix.exp >& $fprefix.out
#  mkdir "PAR_1min_F2030"
#  mv $fprefix.* "PAR_1min_F2030"/.
#  echo "Ended forecast from 20:30 to 21:30 UTC"
#  date

#-----------------------------------------------------------------------------------------

   echo "Running forecast from 20:20 UTC"
   cp ./PAR_1min_F2030/$fprefix.* .
   python run_fcst.py --prefix $fprefix --start 1800 --stop 5400 --nthreads $nthreads --range 1 $nens
   mkdir "PAR_1min_F2020"
   mv $fprefix.* "PAR_1min_F2020"/.
   echo "Ended forecast from 20:20 to 21:20 UTC"
   date

#-----------------------------------------------------------------------------------------

   echo "Running forecast from 20:10 UTC"
   cp ./PAR_1min_F2030/$fprefix.* .
   python run_fcst.py --prefix $fprefix --start 1200 --stop 4800 --nthreads $nthreads --range 1 $nens
   mkdir "PAR_1min_F2010"
   mv $fprefix.* "PAR_1min_F2010"/.
   echo "Ended forecast from 20:10 to 21:10 UTC"
   date

endif
#--------END PAR EXP----------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#  88D experiment

if ( $run_5min == 1 ) then
   echo "STARTING 5 minute experiment.."
   echo "Initializing run with cold start"

#  python pyinit2km_cold.py

#-----------------------------------------------------------------------------------------

#  echo "Running 5 min cycle...."
#  python cold_5min.py -f $fprefix.exp >& $fprefix.out
#  mkdir "PAR_5min_F2030"
#  mv $fprefix.* "PAR_5min_F2030"/.
#  mv Prior* "PAR_5min_F2030"/.
#  mv Inflation* "PAR_5min_F2030"/.
#  echo "Ended 5 min cycle"

#-----------------------------------------------------------------------------------------

   echo "Running forecast from 20:20 UTC"
   cp ./PAR_5min_F2030/$fprefix.* .
   python run_fcst.py --prefix $fprefix --start 1800 --stop 5400 --nthreads $nthreads --range 1 $nens
   mkdir "PAR_5min_F2020"
   mv $fprefix.* "PAR_5min_F2020"/.
   echo "Ended forecast from 20:20 to 21:20 UTC"
   date

#-----------------------------------------------------------------------------------------

   echo "Running forecast from 20:10 UTC"
   cp ./PAR_5min_F2030/$fprefix.* .
   python run_fcst.py --prefix $fprefix --start 1200 --stop 4800 --nthreads $nthreads --range 1 $nens
   mkdir "PAR_5min_F2010"
   mv $fprefix.* "PAR_5min_F2010"/.
   echo "Ended forecast from 20:10 to 21:10 UTC"
   date

endif
#-----END 88D EXP-------------------------------------------------------------------------
