#!/bin/bash
#SBATCH --job-name=dpgen_init
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --partition=hfacnormal01

# load the environment
module purge
source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=/public/home/c1337375425/soft/vasp.5.4.4/bin/:${PATH}
ulimit -s unlimited

# DESC: 此脚本用于批量执行基于ASE DATABASE中的结构计算任务。为一系列ID（从start_id到end_id）执行相应的结构计算流程。每个ID对应一次独立的VASP计算，该计算包括导出POSCAR文件、准备输入文件、运行VASP以及记录日志。需要在当前目录下提供INCAR, POTCAR, 其中K点在INCAR中由KSPACING提供，并且需要配套的python文件用于导出ase.db模块的结构。
# USAGE: sbatch dpgen_init.sh <start_id> <end_id>

dp_export_script="export_fp_stru.py"
log_name="out.log"
start_id=$1
end_id=$2

echo `date` > ${log_name}
echo "author: gchen" >> ${log_name}

echo "" >> ${log_name}
echo "+---------------------------------+" >> ${log_name}
echo "将要计算的id: start: ${start_id}, end: ${end_id} (include end id)" >> ${log_name}
echo "+---------------------------------+" >> ${log_name}
echo "" >> ${log_name}

# 检查是否提供了两个参数
if [ "$#" -ne 2 ]; then
    echo "Error: 需要提供两个参数!" >> ${log_name}
    echo "Usage: $0 <start_id> <end_id>" >> ${log_name}
    exit 1
fi

# 检查id顺序：start_id < end_id
if [ "${start_id}" -gt "${end_id}" ]; then
    echo "Error: 起始值 (${start_id}) 必须小于或等于结束值 (${end_id})!" >> ${log_name}
    exit 1
fi

for (( num=start_id; num<=end_id; num++ )); do
    # echo "当前数字: $num"
    echo "+---------------------------------+" >> ${log_name}
    echo "开始进行 id = ${num} 的结构计算" >> ${log_name}
	echo "Time: $(date +"%Y-%m-%d %H:%M:%S")" >> ${log_name}
    echo "+---------------------------------+" >> ${log_name}
    python export_fp_stru.py ${num}
    if [ $? -ne 0 ]; then
        echo "Error: 导出POSCAR失败!" >> ${log_name}
        echo "" >> ${log_name}
        continue
    fi

    mkdir ${num}

    cp POTCAR INCAR ./${num}
    mv POSCAR ./${num}

    cd ${num}

    srun vasp_std
    if [ $? -eq 0 ]; then
        echo "id = ${num} 的结构计算完成!" >> ../${log_name}
        echo "" >> ../${log_name}
    else
        echo "Error: id = ${num} 的结构计算失败!" >> ../${log_name}
        echo "" >> ../${log_name}
    fi

    cd ..
done