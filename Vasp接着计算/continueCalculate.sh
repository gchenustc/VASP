#!/bin/env bash
#第一个参数：备份后缀，第二个参数：提交任务的命令，第三个参数是否备份OUTCAR，是则填写任意内容
#第二和第三参数可以不填，则只进行备份
subCmd="sbatch"

if ! test $1;then
    $1=1
fi

suffix=$1
script=$2

cp XDATCAR XDATCAR_${suffix}
cp OSZICAR OSZICAR_${suffix}
cp POSCAR POSCAR_${suffix}
cp INCAR INCAR_${suffix}
if test $3;then
    cp OUTCAR OUTCAR_${suffix}
fi

mv CONTCAR POSCAR
if test ${script}; then
    $subCmd ${script}
fi
