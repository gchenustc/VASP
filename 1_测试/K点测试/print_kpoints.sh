#!/bin/env bash

echo "# KPOINTS_VALUE vs TOTAL_ENERGY" >  out.dat
echo "# DATE: $(date)">> out.dat
echo "# PWD": $(pwd) >> out.dat
printf "KPOINTS       Energy\n" >> out.dat


list=("1 1 1" "2 2 2" "3 3 3" "4 4 4" "5 5 5" "6 6 6" "7 7 7" "8 8 8" "9 9 9" "10 10 10" "11 11 11" "12 12 12" "13 13 13")
KpCount=${#list[*]}  # 统计测试的k点个数
for i in $(seq 0 1 $((${KpCount}-1)))
do
    workDir=kpoints_`echo $(echo ${list[$i]}) | awk 'BEGIN{OFS="-"};{print $1,$2,$3}'`
    kpoints=`echo $(echo ${list[$i]}) | awk 'BEGIN{OFS="-"};{print $1,$2,$3}'`
    [ -d $workDir ] || echo "the target dir not existexit"$workDir > error
    energy=`grep -E 'energy[ ]+without[ ]+entropy' $workDir/OUTCAR | tail -1 | awk '{print $4}'`
    printf "%-10s" ${kpoints}"  " >> out.dat
    printf "%15.10f" $energy >> out.dat
    printf "\n" >>out.dat
done
