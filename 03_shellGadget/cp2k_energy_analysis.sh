#!/bin/bash
# 脚本运行时间
echo "Script run time: $(date)" > energy.log

# 作者
echo "Author: gchen" >> energy.log

# 描述和用法
echo "This script finds and extracts the total energy value from .out files in specified subdirectories and plots energy change with primary folder name. Noted that, the user should specify the subdirectory in no.1 parameters. For example: 'bash **.sh calcE'. And then this script will enter the */calcE/*.out and extract the energy into the energy.log" >> energy.log

# 检查参数是否为空
if [ -z "$1" ]; then
    echo "Please provide the name of the subdirectory as an argument."
    exit 1
fi

# 遍历一级文件夹
for dir in */; do
    if [ -d "$dir" ]; then
        sub_dir="$dir$1"
        if [ -d "$sub_dir" ]; then
            file_path=$(find "$sub_dir" -name "*.out")
            echo $file_path
            if [ -n "$file_path" ]; then
                n_atoms=$(grep -o '\- Atoms:[ ].* [0-9].*' "$file_path" | awk '{print $3}')
                energy_value=$(grep -o 'Total energy: [ ].* -[0-9].*' "$file_path" | awk '{print $3}')
                average_energy=$(echo "scale=4; $energy_value/$n_atoms" | bc)
                echo "$sub_dir: $average_energy" >> energy.log
            else
                echo "No .out file found in $sub_dir" >> energy.log
            fi
        fi
    fi
done

# 生成gnuplot数据文件
echo "Generating gnuplot data file..."
awk -F'/' '{print $1, $2, $3}' energy.log > energy_data.txt

# 绘制能量随一级文件夹名称变化的图表
echo "Plotting energy change with primary folder name..."
gnuplot << EOF
set term png
set output 'energy_plot.png'
set xlabel "Primary Folder Name"
set ylabel "Energy Value"
set title "Energy Change with Primary Folder Name"
set xtics rotate by -45
set style data linespoints
set autoscale y
plot "energy_data.txt" using 1:3 with linespoints
EOF
