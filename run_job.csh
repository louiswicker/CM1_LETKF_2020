#!/bin/csh

echo "Starting run_job.csh"

#==> run_Exper: start time for experiment is 2003-05-08 20:40:00

#==> run_Exper: stop  time for experiment is 2003-05-08 21:30:00
#======>>> Step  0:  Run forecast from 2003_05-08 20:40:00 until 2003_05-08 21:00:00
#======>>> Step  1:  Assimilate and run additive noise at 2003_05-08 21:00:00, then run forecast until 2003_05-08 21:05:00
#======>>> Step  2:  Assimilate and run additive noise at 2003_05-08 21:05:00, then run forecast until 2003_05-08 21:10:00
#======>>> Step  3:  Assimilate and run additive noise at 2003_05-08 21:10:00, then run forecast until 2003_05-08 21:15:00
#======>>> Step  4:  Assimilate and run additive noise at 2003_05-08 21:15:00, then run forecast until 2003_05-08 21:20:00
#======>>> Step  5:  Assimilate and run additive noise at 2003_05-08 21:20:00, then run forecast until 2003_05-08 21:25:00
#======>>> Step  6:  Assimilate and run additive noise at 2003_05-08 21:25:00, then run forecast until 2003_05-08 21:30:00
#======>>> Step  7:  Assimilate at 2003_05-08 21:30:00, then run forecast until 2003_05-08 21:45:00

#==> run_Exper: Asynchronous DA is requested, overshoot = 150 sec
#  THIS IS WRONG |


#--------------------------------------------------------------------------

/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 1200 -t 2003,05,08,20,40,00

#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:00:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,00,00 --freq -300

/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,00,00 --write

#--------------------------------------------------------------------------

/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,00,00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,05,00 --freq -300

/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,05,00 --write

#--------------------------------------------------------------------------

/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,05,00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,10,00 --freq -300

/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,10,00 --write

#--------------------------------------------------------------------------

/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,10,00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,15,00 --freq -300

/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,15,00 --write

#--------------------------------------------------------------------------

/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,15,00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,20,00 --freq -300

/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,20,00 --write

#--------------------------------------------------------------------------

/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,20,00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,25,00 --freq -300

/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,25,00 --write

#--------------------------------------------------------------------------

/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,25,00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,30,00 --freq -300
