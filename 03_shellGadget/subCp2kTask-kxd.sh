#!/bin/env bash
# Desc: auto submit the cp2k task, the current folder needs the xx.inp file. Noted, you should to assagin the BASIS,  POTENTIAL and submitted script path, first parameter is .inp file, second is the n.th node, thrid is optional, assagin the basis file, default is BASIS_MOLOPT. For example: ./subCp2kTask-kxd.sh demo.inp 2 BASIS_pob
# Usage: `currentScript.sh xx.inp basis_file_name`, the second parameter could be None, defaule is BASIS_MOLOPT


if [ -n "$3" ];then
    basis=$3
else
    basis=BASIS_MOLOPT
fi

pbs_name_=$1
pbs_name=${pbs_name_%.*}
echo $pbs_name

cp /home/ISSP2/gchenn/Storage/cp2k/$basis ./
cp /home/ISSP2/gchenn/Storage/cp2k/POTENTIAL ./
cp /home/ISSP2/gchenn/Storage/Scripts/submit_scripts/cp2k-node$2.pbs ./

sed -i "s/^#PBS -N .*/#PBS -N $pbs_name/" cp2k-node$2.pbs

qsub cp2k-node$2.pbs
