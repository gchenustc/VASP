#!/bin/env bash
# 优化了文件夹可以为负数
current=-1
for pre in {-20..-60..-20}
do
    next=${pre}

    if [ -d ${next} ];then
        rm -f ./${next}/*
    else
        mkdir ./${next}
    fi
    cp INCAR POSCAR POTCAR ./${next}

    if [ ${current} -ne -1 ];then
        cp ./${current}/CONTCAR ./${next}/POSCAR
    fi

    cd ./${next}
    sed -i /PSTRESS/cPSTRESS=${next}0 INCAR
    srun vasp_std
    current=${next}
    cd ..
done
