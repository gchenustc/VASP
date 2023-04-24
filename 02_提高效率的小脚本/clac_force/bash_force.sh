#!/bin/env bash
# 遍历子文件夹，打印每个文件夹的最大力和平均力，需要getForce.py文件，并指定其路径

echo `date` >> forces.log
echo "energy_density  max_force  ave_force" >> forces.log

for i in `ls`
do
    if [[ -d ${i} ]];then
        cd ${i}
            echo -n "${i}  " >> ../forces.log
            echo `python /public/home/c1337375425/scripts/getForce.py` >> ../forces.log
        cd ..
    fi
done
