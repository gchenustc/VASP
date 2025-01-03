#!/bin/bash
#SBATCH --job-name=hydroStrus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=lzicnormal
#SBATCH --exclude=b01r4n[0208,0302],b01r3n[0208,0304],b01r5n[0204,0303,0402],b01r2n[0408]

# load the environment
module purge
source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=/public/home/c1337375425/soft/vasp.5.4.4/bin/:${PATH}
ulimit -s unlimited

# Author: gchen
# Date: 20250103
# Desc: To generate hydrostatic pressure structures (an accompanied python file is needed, and in the same folder with this script), and calculate the structures. the first position parameter is the scale factor (eg. sbatch *.sh -0.05 means reduce 5% of specified direction, details setting sees the python file). The template POSCAR, POTCAR, INCAR file is needed in the up one level folder. 

# ---------- setting -----------
n=50   # 循环个数
log_name="out.log"
# -------------------------------

echo `date` > ${log_name}
echo "Author: gchen" >> ${log_name}

echo "" >> ${log_name}
echo "+---------------------------------+" >> ${log_name}
echo "循环个数: ${n}; lattice 变化系数: $1" >> ${log_name}
echo "+---------------------------------+" >> ${log_name}
echo "" >> ${log_name}

for ((i=1; i<=n; i++))
do
  echo "+---------------------------------+" >> ${log_name}
  echo "Start no.${i}." >> ${log_name}
  echo "Time: $(date +"%Y-%m-%d %H:%M:%S")" >> ${log_name}
  echo "+---------------------------------+" >> ${log_name}

  # 删除已存在的目录及其内容（如果存在）
  if [ -d "$i" ]; then
    rm -rf "$i"
    echo "Removed and will recreate directory $i" >> ${log_name}
  fi

  # 创建名为$i的目录
  mkdir "$i"
  echo "Now created directory $i." >> ${log_name}

  # 产生结构
  echo "Now create structures." >> ${log_name}
  python gen_hydrostatic_pressure_strus.py $1

  # 创建输入文件
  cp POSCAR ./${i}
  cp ../POTCAR ./${i}
  cp ../INCAR ./${i}

  # 计算
  cd $i
  echo "Start calculation!" >> ../${log_name}
  srun vasp_std
  echo "End calculation!" >> ../${log_name}
  cd ..
done

