# 遍历当前文件夹下的所有文件，生成对应的文件夹，将该文件拷贝进去，再将文件夹attach下的文件拷贝到生成的文件夹下，运行每个文件夹下的vasp_job.sh

declare -A dirs       # 储存POSCAR的名字 dirs[1]=POSCAR-001, dirs[2]=POSCAR-002
no_dir=0               # POSCAR-*的数量，也就是文件夹的数量
for i in `ls`
do
    result=$(echo ${i} | grep "3x1x2")
    if [ "${result}" != "" ];then
        ((no_dir+=1))
        dirs[${i:16:10}]=${i}
    fi
done

for i in ${!dirs[@]}
do
    mkdir ${i}
    cp ${dirs[$i]} ${i}/POSCAR
    cp attach/* ${i}
    cd ${i}
    sbatch vasp_job.sh
    cd ..
done
