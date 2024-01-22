#!/bin/env bash
#Desc: Backup Vasp input params if the task break abruptly. the first param is the suffix for the backup file. Default backup INCAR, POSCAR, XDATCAR, OSZICAR. OUTCAR, if there are other file to backup, write them in the no.2 and later params. eg: ./script.sh 1 POTCAR CHGCAR
#Author: gchen
#Time: 2022

# if not pass any params, set the $1=1
if [ $# -eq 0 ];then
    $1=1
fi

# the no.2 the later params are the extra file wanted to backup.
args=("$@")
for ((i=1; i<${#args[@]}; i++));do
    cp ${args[i]} ${args[i]}_$1
done

cp XDATCAR XDATCAR_$1
cp OSZICAR OSZICAR_$1
cp POSCAR POSCAR_$1
cp INCAR INCAR_$1
cp OUTCAR OUTCAR_$1

cp CONTCAR POSCAR
