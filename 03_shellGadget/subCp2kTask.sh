#!/bin/env bash
# Desc: auto submit the cp2k task, the current folder needs the xx.inp file. Noted, you should to assagin the BASIS,  POTENTIAL and submitted script path
# Usage: `currentScript.sh xx.inp basis_file_name`, the second parameter could be None, defaule is BASIS_MOLOPT


if [ -n "$2" ];then
    basis=$2
else
    basis=BASIS_MOLOPT
fi

cp ~/storage/cp2k/$basis ./
cp ~/storage/cp2k/POTENTIAL ./
cp ~/scripts/sub_scripts/cp2k.slurm ./
sbatch cp2k.slurm $1
