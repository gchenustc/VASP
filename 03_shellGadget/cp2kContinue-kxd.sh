#/bin/env bash
# author: gchen
# time: 20240606

# suffix
NUM=0

# 遍历所有以 "run-" 开头并以 ".log" 结尾的文件
# 寻找已经备份的.log文件的最大后缀为NUM
for file in run-*.log; do
    # 使用basename和awk来分割文件名，获取数字部分
    num=$(basename "$file" | awk -F'[-.]' '{print $3}')
    # 检查数字是否存在
    if [[ $num =~ ^[0-9]+$ ]]; then
        if ((num > NUM));then
            NUM=$num
        fi
    fi
done

NUM=$((NUM+1))

# 备份run-*.log 为 run-*-NUM.log
for file in run-*.log; do
    num=$(basename "$file" | awk -F'[-.]' '{print $3}')
    if ! [[ $num =~ ^[0-9]+$ ]]; then
        mv ${file} ${file%.*}-${NUM}.log
    fi
done


# replace .pbs file
for file in cp2k-node*.pbs; do
    base_name=$(basename "$file" .pbs) # do not contain suffix
    out_file="${base_name}-restart.pbs"
    sed 's/\$PBS_JOBNAME\.inp/\$PBS_JOBNAME-1\.restart/g' ${file}  > ${out_file}
done

qsub cp2k-node*-restart.pbs
