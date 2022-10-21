#!/bin/env bash
incar=INCAR
kpoints=KPOINTS  #要把K点那行数据替换成SET_KPOINTS
poscar=POSCAR
potcar=POTCAR

list=("1 1 1" "2 2 2" "3 3 3" "4 4 4" "5 5 5" "6 6 6" "7 7 7" "8 8 8" "9 9 9" "10 10 10" "11 11 11" "12 12 12" "13 13 13")
KpCount=${#list[*]}  # 统计测试的k点个数
for i in $(seq 0 1 $((${KpCount}-1))) #假如k点数为3，就在0..2之间循环 ##{0..$((${KpCount}-1))} ##{0..$[${KpCount}-1]} ##{0..$((${KpCount}-1))} ##{0..`expr ${KpCount} - 1`}
do
    workDir=kpoints_`echo $(echo ${list[$i]}) | awk 'BEGIN{OFS="-"};{print $1,$2,$3}'`
    if [ ! -d $workDir ]; 
    then 
        mkdir $workDir
    else
        rm $workDir/*
    fi
    
    Kpoints=${list[$i]}
    sed s/SET_KPOINTS/"$Kpoints"/g KPOINTS > $workDir/KPOINTS
    cp $potcar $workDir/POTCAR
    cp $poscar $workDir/POSCAR
    cp $incar $workDir/INCAR
done
