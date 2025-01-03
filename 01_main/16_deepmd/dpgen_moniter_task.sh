#!/bin/env bash
# Author: gchen
# Time: 20241121
# Desc: Check the completion status of tasks in the dpgen working directory. Place this script in the remote computation directory of dpgen and execute it directly. if exist $1=debug, debug mode started.


sub_run_file=$(find . -type f -name "*.sub.run")

for i in ${sub_run_file};do
	echo "start count ${i}"

	# 使用grep查找模式，然后通过sort和uniq去重，最后将结果存储到数组中
	#mapfile -t finished_tag < <(grep -o '[a-f0-9]\{40\}_task_tag_finished' "${i}" | sort | uniq)
	#mapfile -t idx_tag < <(grep -o 'task.[0-9]\{3\}.[0-9]\{6\}' "${i}" | sort | uniq)
	# 读取文件并逐行处理
    # 初始化两个数组
    idx_tag=()
    finished_tag=()

    # 读取文件并逐行处理
    while IFS= read -r line; do
        # 匹配任务目录名（如 task.004.000055）
        if [[ $line =~ cd\ ([^ ]+/)?(task\.[0-9]{3}\.[0-9]{6}) ]]; then
            task_dir=${BASH_REMATCH[2]}
            idx_tag+=("$task_dir")
        fi

        # 匹配任务完成标记文件名（如 9cdd0685529666a406276eaa7a3d000ef216af15_task_tag_finished）
        if [[ $line =~ (!\ -f\ )([a-f0-9]{40}_task_tag_finished) ]]; then
            tag=${BASH_REMATCH[2]]}
            finished_tag+=("$tag")
        fi
    done < "$i"

    # 确保两个数组长度一致
    if [ "${#idx_tag[@]}" -ne "${#finished_tag[@]}" ]; then
        echo "Warning: The number of task directories and finished tags do not match in file ${i}."
    fi

	# 遍历数组，打印每个元素
	#for element in "${finished_tag[@]}"; do
	#    echo "$element"
	#done
	n_total=${#finished_tag[@]}
	n_completed=0
	n_count=0
	for j in "${finished_tag[@]}"; do
		#echo ${finished_tag[${n_count}]}
		#echo ${idx_tag[${n_count}]}
		if find ./${idx_tag[$n_count]} -type f -name "$j" | grep -q .; then
			n_completed=$((n_completed+1))
		elif [ "$1" = "debug" ]; then
				echo "$j not existed, its calculation directly is ${idx_tag[${n_count}]}"
		fi
		n_count=$((n_count+1))
	done
	echo ${i}: total: ${n_total}, now complete ${n_completed} tasks.
	echo
done