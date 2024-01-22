#!/bin/bash

DIR=$(pwd)

nderiv_step_size='2_0.01'

export VASP_RAMAN_RUN='srun vasp_std > vasp.log'

# manually split modes
Modes[0]='01_07'
Modes[1]='08_15'
Modes[2]='16_21'

for m in ${Modes[*]}
do
    VASP_RAMAN_PARAMS="${m}_${nderiv_step_size}"
    export VASP_RAMAN_PARAMS=${VASP_RAMAN_PARAMS}
    #
    FOLDER="Modes_${VASP_RAMAN_PARAMS}"
    #
    echo "Running VASP in ${FOLDER}"
    #
    if [ -d "${FOLDER}" ]; then
        cd "${FOLDER}"
        sbatch raman.slurm
        cd ..
    else
        mkdir "${FOLDER}"
        cd "${FOLDER}"
        ln -s ../OUTCAR.phon ./OUTCAR.phon
        ln -s ../POSCAR.phon ./POSCAR.phon
        ln -s ../POTCAR ./POTCAR
        ln -s ../INCAR ./INCAR
        ln -s ../raman.slurm ./raman.slurm
        ln -s ../vasp_raman.py ./
        #ln -s ../KPOINTS ./KPOINTS
        sbatch  raman.slurm
        cd ..
    fi
done
