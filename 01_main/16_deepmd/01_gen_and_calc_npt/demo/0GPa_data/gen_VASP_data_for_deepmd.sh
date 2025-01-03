#!/bin/bash
#SBATCH --job-name=vasp-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=lzicnormal
#SBATCH --exclude=b01r4n[0208,0302],b01r3n[0304]

# Author: gchen
# Desc: It's for generate the datas for deepmd, you need to create three folder, there are vasp_input_coarse (place INCAR, POTCAR, KPOINTS in it), vasp_input_fine (place INCAR, POTCAR in it, KPOINTS write in the INCAR file), and calculated folder (place this script and poscar in it.). Then set the parameter. There are outer_tem_gap, n_outer, and n_inner following.  outer_tem_gap assigns the temperature for increase, n_outer assigns outer circle times. Positional parameter assigns the pressure. For each inner circle, program will start a coarse NPT calculation and then a fine NPT calculation. Then program will start the outer circle, For outer circle, it would increae the temperature. 
# For example, you want to get the data for POSCAR in 20 GPa and for temperature up to 1000 K, you could set outer_tem_gap=200, n_outer=5, n_inner=5, and run this script. `Sbatch ./Script.py 20`.

# load the environment
module purge
source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=/public/home/c1337375425/soft/vasp.5.4.4/bin/:${PATH}
ulimit -s unlimited

# ===== set =====
pressure=$1
outer_tem_gap=200
n_outer=5    # n_outer plus outer_tem_gap equals the maximum temperature 
n_inner=5

log_name="out.log-1"
# ===== set =====

n_outer=$((n_outer + 1))
n_inner=$((n_inner + 1))

echo "+---------------------------------+" >> ${log_name}
echo `date` > ${log_name}
echo "author: gchen" >> ${log_name}
echo "pressure: ${pressure} GPa" >> ${log_name}
echo "+---------------------------------+" >> ${log_name}
echo "" >> ${log_name}

# 初始化
mkdir ./fine_results
mkdir ./coarse_results

path_coarse="" # 初始化
path_fine="" # 初始化

# 进行不同温度下的计算
# 外循环循环温度
for((i=1;i<${n_outer};i++)); do
    tem=$(( i * outer_tem_gap ))
    echo "+---------------------------------+" >> ${log_name}
    echo "current temperature: ${tem}" >> ${log_name}
    echo "+---------------------------------+" >> ${log_name}

    mkdir ./fine_results/${tem}
    mkdir ./coarse_results/${tem}

    # 内循环 - 分别进行一次粗优化和一次精优化
    for((j=1;j<${n_inner};j++)); do
        echo "current iter: ${j}" >> ${log_name}
        path_coarse_old=${path_coarse}
        path_fine_old=${path_fine}
        path_coarse="./coarse_results/${tem}/${j}"
        path_fine="./fine_results/${tem}/${j}"
        mkdir ${path_coarse}
        mkdir ${path_fine}

        # ====== 粗优化 ====== 
        echo "current is coarse calculations." >> ${log_name}
        cp ../vasp_input_coarse/INCAR ${path_coarse}/
        cp ../vasp_input_coarse/POTCAR ${path_coarse}/
        cp ../vasp_input_coarse/KPOINTS ${path_coarse}/
        
        if [ ${i} -eq 1 ] && [ ${j} -eq 1 ];then
            cp ./POSCAR ${path_coarse}/
            echo 1
        else
            cp ${path_fine_old}/CONTCAR ${path_coarse}/POSCAR
            echo 2
        fi

        cd ${path_coarse}

        # 替换INCAR
        # 改变压力
        sed -i 's/PSTRESS.*=.*/PSTRESS = '${pressure}'0/' INCAR
        # 改变温度
        sed -i 's/TEBEG.*=.*/TEBEG = '${tem}'/' INCAR
        sed -i 's/TEEND.*=.*/TEEND = '${tem}'/' INCAR

        # 执行计算**
        #cp POSCAR CONTCAR
        srun vasp_std

        cd - &> /dev/null

        # ====== 精优化 ====== 
        echo "current is fine calculations." >> ${log_name}
        echo "" >> ${log_name}
        cp ../vasp_input_fine/INCAR ${path_fine}/
        cp ../vasp_input_fine/POTCAR ${path_fine}/
        
        cp ${path_coarse}/CONTCAR ${path_fine}/POSCAR

        cd ${path_fine}/

        # 替换INCAR
        # 改变压力
        sed -i 's/PSTRESS.*=.*/PSTRESS = '${pressure}'0/' INCAR
        # 改变温度
        sed -i 's/TEBEG.*=.*/TEBEG = '${tem}'/' INCAR
        sed -i 's/TEEND.*=.*/TEEND = '${tem}'/' INCAR

        # 执行计算**
        #cp POSCAR CONTCAR
        srun vasp_std

        cd - &> /dev/null
        
    done
done


