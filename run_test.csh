#!/bin/csh
echo "Starting runit.csh"

set dir = "test2kmB"

echo "running job "$dir

date

run_exper.py -d RUN_LETKF -i >& out_$dir
mv RUN_LETKF $dir

date

ens.py -d $dir -t 2003,5,8,21,0,0 -v W --plot9 #--zoom 80 160 80 160
ens.py -d $dir -t 2003,5,8,21,20,0 -v W --plot9 #--zoom 80 160 80 160
ens.py -d $dir -t 2003,5,8,21,40,0 -v W --plot9 #--zoom 80 160 80 160
ens.py -d $dir -t 2003,5,8,22,0,0 -v W --plot9 #--zoom 80 160 80 160
#
ens.py -d $dir -t 2003,5,8,21,0,0 -v DBZ --plot8
ens.py -d $dir -t 2003,5,8,21,20,0 -v DBZ --plot8 
ens.py -d $dir -t 2003,5,8,21,40,0 -v DBZ --plot8 # --zoom 80 160 80 160
ens.py -d $dir -t 2003,5,8,22,0,0 -v DBZ --plot8 #--zoom 80 160 80 160

#vort_prob_CM1.py -d $dir -t 2011,5,24,20,35,0 --fcst 3000 60 --noshow
#
#cp $dir*.pdf $dir/Plots/.

python DBZ_CR.py -d $dir -t DBZ_CR_$dir --noshow
python DBZ_INV.py -d $dir -t DBZ_INV_$dir --noshow

python VR_CR.py -d $dir -t VR_CR_$dir --noshow
python VR_INV.py -d $dir -t VR_INV_$dir --noshow

#mv VR_*.pdf $dir/Plots/.

date
exit(0)
#---------------------------------------------------------
