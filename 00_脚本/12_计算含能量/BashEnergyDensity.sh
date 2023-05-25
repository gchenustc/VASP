#!/bin/env bash
key_words=3P  # 遍历子文件夹，当文件夹名含有该关键词则为候选文件夹，在脚本中要给定计算能量密度的文件
# 要在脚本中给定CalcEnergyDensity.py的位置

echo `date` >> energy_density.log
#echo "folder  energy_density  max_force  ave_force" >> energy_density.log
echo "folder  energy_density" >> energy_density.log

for i in `ls`
do
    if [[ -d ${i} ]] && [[ ${i} =~ ${key_words} ]];then
        cd ${i}
            echo -n "${i}  " >> ../energy_density.log
            echo -n "`python /public/home/c1337375425/scripts/CalcEnergyDensity.py`  " >> ../energy_density.log
            #echo `python /public/home/c1337375425/scripts/getForce.py` >> ../energy_density.log
        cd ..
    fi
done
