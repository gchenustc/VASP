#!/bin/env bash
#SBATCH --job-name=Relax_N2F2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=hfacnormal04

# description: 变压逐步优化，需要INCAR,POSCAR,POTCAR,vasp_job.sh文件。注意需要设置KSPACING

# load the environment
module purge
source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=/public/home/c1337375425/soft/vasp.5.4.4/bin/:${PATH}

current=-1
# 从160到0，步长为20遍历
for pre in {160..0..20}
do
    next=${pre}

    # 检查是否存在文件夹${next}，如果存在则清空文件夹内的文件，否则创建文件夹
    if [ -d ${next} ]; then
        rm -f ./${next}/*
    else
        mkdir ./${next}
    fi

    # 复制文件到${next}文件夹
    cp INCAR POSCAR POTCAR ./${next}

    # 如果current不等于-1，则将CONTCAR文件从${current}文件夹复制到${next}文件夹下的POSCAR文件中
    if [ ${current} -ne -1 ]; then
        cp ./${current}/CONTCAR ./${next}/POSCAR
    fi

    # 进入${next}文件夹
    cd ./${next}
    # 在INCAR文件中替换PSTRESS的值为${next}0
    sed -i 's/PSTRESS=.*/PSTRESS='${next}'0/' INCAR
    # 执行vasp_std命令
    srun vasp_std
    # 返回上一级目录
    cd ..

    # 更新current的值为${next}
    current=${next}
done

