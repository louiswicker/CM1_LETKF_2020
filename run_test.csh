#!/bin/csh
echo "Starting run_test.csh"

set dir = "test"

echo "running job "$dir

date

create_run_letkf.py >& out_$dir
run_fcst.py -e RUN_LETKF/RUN_LETKF.exp -i >& out_$dir
ens.py -e RUN_LETKF/RUN_LETKF.exp --init0 -t 2003,5,8,20,40,0 --write >& out_$dir

# run job script

run_job.csh >& out_$dir

cp -R RUN_LETKF $dir

date

ens.py -d RUN_LETKF -t 2003,5,8,21,0,0 -v W --plot9 #--zoom 80 160 80 160
ens.py -d RUN_LETKF -t 2003,5,8,21,15,0 -v W --plot9 #--zoom 80 160 80 160
ens.py -d RUN_LETKF -t 2003,5,8,21,30,0 -v W --plot9 #--zoom 80 160 80 160
ens.py -d RUN_LETKF -t 2003,5,8,22,0,0 -v W --plot9 #--zoom 80 160 80 160
#
ens.py -d RUN_LETKF -t 2003,5,8,21,0,0 -v DBZ --plot8
ens.py -d RUN_LETKF -t 2003,5,8,21,15,0 -v DBZ --plot8 
ens.py -d RUN_LETKF -t 2003,5,8,21,30,0 -v DBZ --plot8 # --zoom 80 160 80 160
ens.py -d RUN_LETKF -t 2003,5,8,22,0,0 -v DBZ --plot8 #--zoom 80 160 80 160

#vort_prob_CM1.py -d $dir -t 2011,5,24,20,35,0 --fcst 3000 60 --noshow
#
#cp $dir*.pdf $dir/Plots/.

python DBZ_CR.py -d RUN_LETKF -t DBZ_CR_$dir --noshow
python DBZ_INV.py -d RUN_LETKF -t DBZ_INV_$dir --noshow

python VR_CR.py -d RUN_LETKF -t VR_CR_$dir --noshow
python VR_INV.py -d RUN_LETKF -t VR_INV_$dir --noshow

#mv VR_*.pdf $dir/Plots/.

date
exit(0)
#---------------------------------------------------------
