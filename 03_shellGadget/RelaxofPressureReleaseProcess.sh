#!/bin/env bash
#SBATCH --job-name=relax
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=hfacnormal04

module purge
source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=/public/home/c1337375425/soft/vasp.5.4.4/bin/:${PATH}

################################################################################
# Description:
# This script performs a step-by-step relaxation process using VASP. It iteratively
# optimizes the structures from an initial configuration to a final configuration
# by gradually reducing(increasing) the pressure.
#
# Usage:
# 1. Place this script in the same directory as the VASP input files (INCAR(KSPACING setting is necessary), POSCAR, POTCAR).
# 2. Modify the parameters (e.g., initial pressure, pressure step size) if needed.
# 3. Run the script using the command: bash vasp_relax.sh
# 4. If the task is interrupted, parameter "current" is set to the folder where the last step has been calculated and meanwhile change the initial_pressure and pressure_step. For example, "initial_pressure=160", "pressure_step=20", but task interrupted for 100, so I set the "current=120", "initial_pressure=100", and "pressure_step=20".
################################################################################

###### Set the initial pressure and pressure step size ######
initial_pressure=160
pressure_step=20
###### ----------------------------------------------- ######

# current folder equals -1 means there are no any calculations now.
# If the task is interrupted, current is set to the folder where the last step "has been calculated"
current=-1
for pre in $(seq ${initial_pressure} -${pressure_step} 0)
do
    next=${pre}

    if [ -d ${next} ]; then
        rm -f ./${next}/*
    else
        mkdir ./${next}
    fi

    cp INCAR POSCAR POTCAR ./${next}

    # 如果current不等于-1，则将CONTCAR文件从${current}文件夹复制到${next}文件夹下的POSCAR文件中
    if [ ${current} -ne -1 ]; then
        cp ./${current}/CONTCAR ./${next}/POSCAR
    fi

    cd ./${next}

    # Replace the PSTRESS value in the INCAR file with the current pressure
    sed -i 's/PSTRESS=.*/PSTRESS='${next}'0/' INCAR
    
    srun vasp_std
    
    cd ..

    # Update the current directory
    current=${next}
done

