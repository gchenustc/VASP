#!/bin/bash
#SBATCH --job-name=NF1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=hfacnormal02

# 用途：不同压力下的结构优化，文件夹名是 kbar 值

# load the environment
module purge
source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=/public/home/c1337375425/soft/vasp.5.4.4/bin/:${PATH}


# 每一个循环提供三个参数，比如"2000 7 1E-5"指的是 pstress=2000，isif=3，ediffg=1E-5
arrays=("2000 7 1E-5" "1800 7 1E-5" "1600 7 1E-5" "1400 7 1E-5" "1200 7 1E-5" "1000 3 -1E-3" "800 3 -1E-3" "600 3 -1E-3" "500 3 -1E-3" "400 3 -1E-3" "300 3 -1E-3" "250 3 -1E-3" "200 3 -1E-3" "150 3 -1E-3" "100 3 -1E-3" "50 3 -1E-3" "0 3 -1E-3")


n_dirs=$((${#arrays[*]}-1)) # 0 is begin

for n in `seq 0 1 ${n_dirs}`
do
    info_array=(${arrays[$n]})
    pstress=${info_array[0]}
    isif=${info_array[1]}
    ediffg=${info_array[2]}

    mkdir ${pstress}
    if [ ${n} -eq 0 ]; then
        cp INCAR POTCAR POSCAR ${pstress}
    else
        info_array_last=(${arrays[${n}-1]})
        pstress_last=${info_array_last[0]}
        cp INCAR POTCAR ${pstress}
        cp ${pstress_last}/CONTCAR ${pstress}/POSCAR
    fi

    cd ${pstress}
	
	# 更改压力
    sed -i -r s/pstress.*/pstress="${pstress}"/1 INCAR
	# 更改isif
    sed -i -r s/isif.*/isif="${isif}"/1 INCAR
	# 更改ediffg
    sed -i -r s/ediffg.*/ediffg="${ediffg}"/1 INCAR
    srun vasp_std

    cd ..
done
