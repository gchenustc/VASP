#!/bin/env bash

echo "# ENCUT_VALUE vs TOTAL_ENERGY" >  out.dat
echo "# DATE: $(date)">> out.dat
echo "# PWD": $(pwd) >> out.dat

printf " ENCUT         Energy\n" >> out.dat
####################需要修改的参数###########################
list="300 320 340 360 380 400 420 440 460 480 500 520 540 560 580 600 620 640 660 680 700 720 740 760 780 800"
####################需要修改的参数###########################
for i in $list
do
    workDir=Encut_$i
    [ -d $workDir ] || echo "the target dir not existexit"\>$workDir > error
    energy=`grep -E 'energy[ ]+without[ ]+entropy' $workDir/OUTCAR | tail -1 | awk '{print $4}'`
    printf "%-7.2f" $i >> out.dat
    printf "%20.10f" $energy >> out.dat
    printf "\n" >>out.dat
done
