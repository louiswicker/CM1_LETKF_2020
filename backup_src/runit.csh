#!/bin/csh
create_run_letkf.py
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp -i
python ens.py -e RUN_LETKF/RUN_LETKF.exp --init0 -t 2003,5,8,20,40,0 --write
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 1800 --time 2003,5,8,20,40,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,10,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,10,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,15,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,15,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,20,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,20,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,25,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,25,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,30,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,30,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,35,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,35,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,40,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,40,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,45,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,45,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,50,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,50,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,21,55,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,21,55,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,22,0,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 300 --time 2003,5,8,22,0,0

python run_filter --exper RUN_LETKF/RUN_LETKF.exp --time 2003,05,8,22,5,0 
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp --run_time 1800 --time 2003,5,8,22,05,0

