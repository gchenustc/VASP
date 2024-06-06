# !/bin/env bash
#author: gchen
#time: 20240517
#desc: 将cp2k结构优化（只适用于不优化胞的结构优化）或者分子动力学产生的.xyz文件转换为.vasp文件，其中脚本的参数：$1 = *.inp，$2 = *.xyz，如果存在$3，则该值指的是步数，即仅仅产生用户指定步数的结构。


# 胞的基矢
cell_str=$(grep -E '&CELL|&cell' -A 5 input.inp | grep -E '^\s*[A|B|C]\s+-?[0-9]+\.[0-9]+[Ee]?[+-]?[0-9]?+\s+-?[0-9]+\.[0-9]+[Ee]?[+-]?[0-9]?+\s+-?[0-9]+\.[0-9]+[Ee]?[+-]?[0-9]?+$' | awk '{print $2, $3, $4}')

read -r -a cell_arr <<< ${cell_str}

# 每一步开始的行号
ln=$(awk '{ 
    if ($0 ~ /i = *[0-9]+, E = *-?[0-9]+\.[0-9]+/) {
        print NR
    }
}' $2)

read -r -a ln_arr <<< ${ln}

#echo "${#ln_arr[@]}"

n_atoms=$(sed -n '1p' $2)

# echo "$n_atoms"

# 储存步数，能量的信息，索引是步数，值是在.xyz中开始的行号
declare -A main_arr

for ln in ${ln_arr[@]}
do
    ln_info=$(sed -n "${ln}p" $2)
    i_value=$(echo $ln_info | awk '{print $3}' | sed 's/.$//') 
    E_value=$(echo $ln_info | awk '{print $6}')
    main_arr["$i_value"]=$ln
done

parent_folder=each_step_structures
if [ -e ${parent_folder} ];then
    rm -rf ${parent_folder}
fi
mkdir ${parent_folder}
cd ${parent_folder}

idx=0
for key in ${!main_arr[@]}
do  
    #echo "${main_arr[$key]}"

    if [ -n "$3" ];then
        if [ $key != "$3" ];then
            continue
        fi
    fi

    ln=${main_arr[$key]}
    start_ln=$((ln+1))
    end_ln=$((ln+n_atoms))
    #echo "$start_ln, $end_ln"
    coord_arr_temp=($(sed -n "${start_ln},${end_ln}p" ../$2 | sort))
    #echo ${coord_arr_temp[0]}

    coord_arr=()  # 每一行的坐标
    declare -A count_n_element  # key是元素名，value是该元素的数量
    
    idx_=0
    coord_each_line=""  # 用来加和坐标，如果坐标满一行，则情况重新加和


    # 用来统计元素数量
    total=${#coord_arr_temp[@]}
    for((i=0;i<$total;i++)); do

        coord=${coord_arr_temp[$i]}

        # 如果 ! i%4，则coord是元素标签
        if ! (( i % 4 )); then
            coord_arr[$idx_]=$coord_each_line
            coord_each_line=""

            if [ -n ${count_n_element[$coord]} ];then
                count_n_element[$coord]=$(( count_n_element[$coord]+1 ))
            else
                count_n_element[$coord]=0
            fi

            idx_=$((idx_+1))

        else
            coord_each_line="${coord_each_line}  $coord"
        fi

        if (( i == total-1 )); then
            coord_arr[$idx_]=$coord_each_line
            coord_each_line=""
        fi
    done

    #for i in ${count_n_element[@]}; do
    #    echo "$i"
    #done


    echo "generate by gchen" > ${key}.vasp
    echo "   1.00000000" >> ${key}.vasp
    # 晶格常数
    temp_i=0
    for i in ${cell_arr[@]};do
        if ! (( temp_i % 3 )) && ((temp_i != 0)) ;then
            echo "" >> ${key}.vasp
        fi
        echo -n "      ${i}" >> ${key}.vasp
        temp_i=$((temp_i+1))
    done

    # 元素
    echo "" >> ${key}.vasp
    for k in ${!count_n_element[@]};do
        echo -n "   $k" >> ${key}.vasp
    done

    # 元素数量
    echo "" >> ${key}.vasp
    for value in ${count_n_element[@]};do
        echo -n "   $value" >> ${key}.vasp
    done

    echo "" >> ${key}.vasp
    echo "C"  >> ${key}.vasp

    # 原子坐标
    temp_i=0
    for value in ${coord_arr[*]};do
        if ! (( temp_i % 3 )) && ((temp_i != 0)) ;then
            echo "" >> ${key}.vasp
        fi
        echo -n "   $value" >> ${key}.vasp
        temp_i=$((temp_i+1))
    done


    idx=$((idx+1))

    # 清空
    for key in ${!count_n_element[@]};do
        unset count_n_element[$key]
    done

done



