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

#==> run_Exper: Now   Model TIME is  0  seconds 

#==> run_Exper: Later Model TIME is  1200.0  seconds 

#==> run_Exper: Step  0:  Run forecast from 2003_05-08 20:40:00 until 2003_05-08 21:00:00

#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_20:40:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 1200 -t 2003,05,08,20,40,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:00:00


#--------------------------------------------------------------------------

#==> run_Exper: Now   Model TIME is  1200  seconds 

#==> run_Exper: Later Model TIME is  1650.0  seconds 

#==> run_Exper: Step  1:  Assimilate and run additive noise at 2003_05-08 21:00:00, then run forecast until 2003_05-08 21:05:00

#==> run_Exper:  Assimilation is being called, the time:  2003-05-08 21:00:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,00,00 --freq -300


#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:00:00

#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:00:00

#==> run_Exper:  Additive noise is based on composite reflectivity


/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,00,00 --write


#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_21:00:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,00,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:05:00


#--------------------------------------------------------------------------

#==> run_Exper: Now   Model TIME is  1500  seconds 

#==> run_Exper: Later Model TIME is  1950.0  seconds 

#==> run_Exper: Step  2:  Assimilate and run additive noise at 2003_05-08 21:05:00, then run forecast until 2003_05-08 21:10:00

#==> run_Exper:  Assimilation is being called, the time:  2003-05-08 21:05:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,05,00 --freq -300


#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:05:00

#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:05:00

#==> run_Exper:  Additive noise is based on composite reflectivity


/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,05,00 --write


#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_21:05:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,05,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:10:00


#--------------------------------------------------------------------------

#==> run_Exper: Now   Model TIME is  1800  seconds 

#==> run_Exper: Later Model TIME is  2250.0  seconds 

#==> run_Exper: Step  3:  Assimilate and run additive noise at 2003_05-08 21:10:00, then run forecast until 2003_05-08 21:15:00

#==> run_Exper:  Assimilation is being called, the time:  2003-05-08 21:10:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,10,00 --freq -300


#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:10:00

#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:10:00

#==> run_Exper:  Additive noise is based on composite reflectivity


/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,10,00 --write


#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_21:10:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,10,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:15:00


#--------------------------------------------------------------------------

#==> run_Exper: Now   Model TIME is  2100  seconds 

#==> run_Exper: Later Model TIME is  2550.0  seconds 

#==> run_Exper: Step  4:  Assimilate and run additive noise at 2003_05-08 21:15:00, then run forecast until 2003_05-08 21:20:00

#==> run_Exper:  Assimilation is being called, the time:  2003-05-08 21:15:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,15,00 --freq -300


#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:15:00

#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:15:00

#==> run_Exper:  Additive noise is based on composite reflectivity


/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,15,00 --write


#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_21:15:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,15,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:20:00


#--------------------------------------------------------------------------

#==> run_Exper: Now   Model TIME is  2400  seconds 

#==> run_Exper: Later Model TIME is  2850.0  seconds 

#==> run_Exper: Step  5:  Assimilate and run additive noise at 2003_05-08 21:20:00, then run forecast until 2003_05-08 21:25:00

#==> run_Exper:  Assimilation is being called, the time:  2003-05-08 21:20:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,20,00 --freq -300


#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:20:00

#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:20:00

#==> run_Exper:  Additive noise is based on composite reflectivity


/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,20,00 --write


#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_21:20:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,20,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:25:00


#--------------------------------------------------------------------------

#==> run_Exper: Now   Model TIME is  2700  seconds 

#==> run_Exper: Later Model TIME is  3150.0  seconds 

#==> run_Exper: Step  6:  Assimilate and run additive noise at 2003_05-08 21:25:00, then run forecast until 2003_05-08 21:30:00

#==> run_Exper:  Assimilation is being called, the time:  2003-05-08 21:25:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,25,00 --freq -300


#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:25:00

#==> run_Exper:  Additive noise is being called, the time: 2003-05-08 21:25:00

#==> run_Exper:  Additive noise is based on composite reflectivity


/Users/Louis.Wicker/cm1_letkf_2020/ens.py -e RUN_LETKF/RUN_LETKF.exp --crefperts -t 2003,05,08,21,25,00 --write


#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_21:25:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 -t 2003,05,08,21,25,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:30:00


#--------------------------------------------------------------------------

#==> run_Exper: Now   Model TIME is  3000  seconds 

#==> run_Exper: Later Model TIME is  3900  seconds 

#==> run_Exper: Step  7:  Assimilate at 2003_05-08 21:30:00, then run forecast until 2003_05-08 21:45:00

#==> run_Exper:  Assimilation is being called, the time:  2003-05-08 21:30:00

/Users/Louis.Wicker/cm1_letkf_2020/run_filter.py --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,08,21,30,00 --freq -300


#==> run_Exper: Integrating ensemble members from date and time: 2003-05-08_21:30:00


/Users/Louis.Wicker/cm1_letkf_2020/run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 1800 -t 2003,05,08,21,30,00


#==> run_Exper: Integrated ensemble members to time: 2003-05-08_21:45:00


#==> run_Exper:  Elapsed WALL CLOCK TIME FOR TOTAL EXPERIMENT: 0.001675 


#==> run_Exper:  Elapsed WALL CLOCK TIME FOR FILTER:           0.000385 


#==> run_Exper:  Elapsed WALL CLOCK TIME FOR MODEL RUNS:       0.000593 


#==> run_Exper:  Elapsed WALL CLOCK TIME FOR ADDITIVE NOISE:   0.000501 

