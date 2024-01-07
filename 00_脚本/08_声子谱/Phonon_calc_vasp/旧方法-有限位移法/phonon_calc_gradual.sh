#!/bin/bash
#SBATCH --job-name=phonon
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=hfacnormal03

# load the environment
module purge
source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=/public/home/c1337375425/soft/vasp.5.4.4/bin/:${PATH}

#描述：本脚本是利用phonopy结合vasp软件，用有限位移法计算phonon的全自动处理程序，提交任务的作业管理系统是sbatch系统，一共只提交一个命令
#需要提供的文件：POSCAR(未扩胞的)，INCAR(要有KSPACING和KGAMMA参数)，POTCAR，作业管理系统提交vasp命令的脚本

#**********
phonopy -d --dim="3 1 3" -c POSCAR  # 设定超胞倍数
#**********

declare -A dirs       # 储存POSCAR的名字 dirs[1]=POSCAR-001, dirs[2]=POSCAR-002
n_dir=0               # POSCAR-*的数量，也就是文件夹的数量
for i in `ls`
do
    result=$(echo ${i} | grep "POSCAR-")
    if [ "${result}" != "" ];then
        ((n_dir+=1))
        dirs[$n_dir]=${i}
    fi
done

for i in `seq 1 1 ${n_dir}`
do
mkdir disp-$i
cp ${dirs[${i}]} disp-${i}/POSCAR
cp INCAR disp-$i/INCAR
cp POTCAR disp-$i/POTCAR
cp ${sub_script} disp-$i/sub.sh

cd disp-$i
    srun vasp_std
cd ..
done 
