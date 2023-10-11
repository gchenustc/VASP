#!/bin/env bash
# 在results里用cak.py --vasp命令，产生一堆dir_**文件
# 将这些文件拷贝到另一个文件夹，脚本在此文件夹运行，在该文件夹下还要准备INCAR，POTCAR和submit.sh

for i in *;do
    if [ -d ${i} ] && [[ ${i} = dir_* ]];then

        cd ${i}
        cd dir_0.1

        for ((j=1; j<=5; j++)); do
            if [ -f UCell_${j}_* ];then
                mkdir ${j}
                cp ../../INCAR ${j}
                cp ../../POTCAR ${j}
                cp ../../submit.sh ${j}
                cp UCell_${j}_* ${j}/POSCAR
                cd ${j}
                sbatch submit.sh
                cd ..
            fi
        done

        cd ..
        cd ..


    fi
done
