#!/bin/env bash
incar=INCAR
kpoints=KPOINTS  #要把K点那行数据替换成SET_KPOINTS
poscar=POSCAR
potcar=POTCAR

######################需要修改的参数########################
list="300 320 340 360 380 400 420 440 460 480 500 520 540 560 580 600 620 640 660 680 700 720 740 760 780 800"
######################需要修改的参数########################
for i in $list
do
    workDir=Encut_$i
    if [ ! -d $workDir ];
    then
        mkdir $workDir
    else
        rm $workDir/*
    fi

    sed s/SET_ENCUT/"ENCUT = $i"/g INCAR > $workDir/INCAR
    cp $potcar $workDir/POTCAR
    cp $poscar $workDir/POSCAR
    cp $kpoints $workDir/KPOINTS
done
