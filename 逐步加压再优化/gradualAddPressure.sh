#!/bin/bash
#SBATCH -p v5_192
#SBATCH -N 1
#SBATCH -n 48
source /public1/soft/modules/module.sh
module load  mpi/intel/17.0.5-cjj
export PATH=/public1/home/scfa1114/soft:$PATH

pressure="20 40 60 80 100"
current=0
for pre in ${pressure}
do
    next=${pre}

    if [ -d ${next} ];then
        rm -f ${next}/*
    else
        mkdir ${next}
    fi
    cp INCAR POSCAR POTCAR KPOINTS sub.sh ./${next}

    if [ ${current} -ne 0 ];then
        cp ${current}/CONTCAR ${next}/POSCAR
    fi

    cd ${next}
    sed -i /PSTRESS/cPSTRESS=${next}0 INCAR
    srun vasp_std
    current=${next}
    cd ..
done

