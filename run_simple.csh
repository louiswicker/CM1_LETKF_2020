#!/bin/csh

echo "Starting run_simple_bench.csh"

set dir = "simple_exper"

echo "running job "$dir

date

# run job script

setenv PYTHONUNBUFFERED TRUE

python run_simple_exper.py >& $dir.out

cp -R RUN_LETKF $dir

cp *.pdf $dir/Plots/.

date
exit(0)
#---------------------------------------------------------
