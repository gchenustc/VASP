#!/bin/sh
#描述：本脚本是利用phonopy结合vasp软件，用有限位移法计算phonon的全自动处理程序，提交任务的作业管理系统是sbatch系统，有限位移法会创建多个文件夹，有多少个文件夹就同时提交多少个任务。
#需要提供的文件：POSCAR(未扩胞的)，INCAR(要有KSPACING和KGAMMA参数)，POTCAR，作业管理系统提交vasp命令的脚本

#**********
sub_script=vasp_job.sh  # 提交任务的脚本名
phonopy -d --dim="3 1 3" -c POSCAR  # 设定超胞倍数
wait
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
    # sed -i "s/#SBATCH --job-name=.*/&${i}/1" sub.sh   # 用序号作为作业名
    sbatch sub.sh
#wait
cd ..
done 
