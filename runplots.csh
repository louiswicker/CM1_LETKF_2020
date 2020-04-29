#!/bin/csh
#
set dir = "PAR_A5_THP"

date

ens.py -d $dir -t 2011,5,24,20,40,0 -v W --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,20,50,0 -v W --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,0,0 -v W --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,10,0 -v W --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,20,0 -v W --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,30,0 -v W --plot9 --zoom 80 160 80 160

ens.py -d $dir -t 2011,5,24,20,40,0 -v DBZ --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,20,50,0 -v DBZ --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,0,0 -v DBZ --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,10,0 -v DBZ --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,20,0 -v DBZ --plot9 --zoom 80 160 80 160
ens.py -d $dir -t 2011,5,24,21,30,0 -v DBZ --plot9 --zoom 80 160 80 160

cp $dir*.pdf $dir/Plots/.

#python VR_CR.py -d $dir -t VR_CR_$dir
#python VR_INV.py -d $dir -t VR_INV_$dir
#
#cp VR_*.pdf $dir/Plots/.

date

