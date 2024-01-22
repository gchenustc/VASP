#!/bin/env bash
#Desc: After completing the tasks of calculating phonon, then use this script. the script accepts params which are the names of folders, the programmers will enter appointed folders and execute the process of post-processing. Need phonopy and vaspkit package
#Author: gchen
#Time: 2023-11-21 19:34

for i in "$@"
do
    path="${PWD}/${i}"
    cd ${path}

    # 调用vaspkit得到能带高对称点
    echo -e "305\n3" | vaspkit

    # 编辑KPATH.in文件
    mv KPATH.phonopy band.conf
    sed -i "8s/^/#/" band.conf
    sed -i "2s/DIM =/DIM = 1 1 1/" band.conf
    sed -i "1s/.*/NPOINTS = 100/" band.conf
    sed -i "6s/.*/MP = 10 10 10/" band.conf 

    # 调用phonopy画图
    phonopy --fc vasprun.xml
    phonopy -c POSCAR band.conf -p -s
    phonopy-bandplot --gnuplot > phonon.out

    cd ..
done
