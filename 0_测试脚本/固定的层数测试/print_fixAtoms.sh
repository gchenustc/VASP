#!/bin/bash
echo "# Fixed_atoms_number vs TOTAL_ENERGY" >  out.dat
echo "# DATE: $(date)">> out.dat
echo "# PWD": $(pwd) >> out.dat
printf "Fixed_number       Energy\n" >> out.dat


for fix_num in $fix_list
do
    # 创建文件夹
    workDir=fix-${fix_num}atoms
    [ -d $workDir ] || echo "the target dir not existexit"-->$workDir
    energy=`grep -E 'energy[ ]+without[ ]+entropy' $workDir/OUTCAR | tail -1 | awk '{print $4}'`
    printf "%-10s" ${fix_num}"  " >> out.dat
    printf "%15.10f" $energy >> out.dat
    printf "\n" >>out.dat
done
