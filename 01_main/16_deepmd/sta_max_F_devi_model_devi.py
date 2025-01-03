import os
import glob
import argparse
import json
import numpy as np
"""
DESC: This script is used to calculate the maximum force deviation of the structures generated in the model_deviation step of dpgen. Place this script in the 01.model_devi directory and execute it directly.
"""

path_pattern="task.00*/job*"
# 读取 温度
with open('cur_job.json', 'r') as file:
    data = json.load(file)
    temps = data['temps']
    print(f"包含的温度: {temps}")
# temps = [100,200,300,400,500]

def is_number(s):
    try:
        float(s)  # 尝试将 s 转换为浮点数
        return True
    except ValueError:
        return False

def paste_files(file_list, output_file):
    # 用来存储每个文件的内容
    columns = []

    # 读取每个文件的内容
    for file_name in file_list:
        with open(file_name, 'r') as f:
            # 将文件内容按行读取并移除空白字符
            column = [line.strip() for line in f.readlines()]
            columns.append(column)

    # 确定结果的行数 (所有文件应该行数相同)
    num_rows = max(len(column) for column in columns)

    # 创建一个列表来存储输出的行
    output_lines = []

    for i in range(num_rows):
        # 组合每列的内容，处理可能不同长度的问题
        row = []
        for column in columns:
            # 如果当前列有此行，添加内容；否则添加空字符串
            row.append(column[i] if i < len(column) else '')

        # 用制表符连接每列
        output_lines.append('\t'.join(row))

    # 写入输出文件
    with open(output_file, 'w') as out_f:
        out_f.write('\n'.join(output_lines))

for t in temps:
    # save dirs_in_this_temps:
    path_dirs = "dirs-"+ str(t)+"K"
    if os.path.exists(path_dirs):
        os.remove(path_dirs)

    # 查找包含特定温度的行，并写入 fdr.dat
    for file_name in glob.glob(path_pattern):
        with open(file_name, 'r') as file:
            for line in file:
                if f'"temps": {t}' in line:
                    print(f"当前文件夹路径: {file_name}，该路径的lammps 在 {t}K 下计算.")
                    dir_name = os.path.dirname(file_name)

                    # 写入 fdr.dat，去除重复项

                    # 首先以读取模式打开文件
                    already_exists = False
                    if os.path.exists(path_dirs):
                        with open(path_dirs, 'r') as fdr_file:
                            content = fdr_file.read()
                            if dir_name in content:
                                already_exists = True

                    # 如果目录名不存在，则写入文件
                    if not already_exists:
                        with open(path_dirs, 'a') as fdr_file:
                            fdr_file.write(dir_name + '\n')

# 删除 ${t}K 文件
    tK_file = f"{t}K"
    if os.path.exists(tK_file):
        os.remove(tK_file)

# 写入legend
    with open(tK_file, 'a') as temp_file:
        temp_file.write(str(t)+ 'K' + '\n')

# 收集数据
    with open(path_dirs, 'r') as fdr_file:
        for line in fdr_file:  # eg: task.000.000000
            dir_name = line.strip()
            model_devi_file = os.path.join(dir_name, 'model_devi.out')
            if os.path.exists(model_devi_file):
                with open(model_devi_file, 'r') as model_file:
                    for model_line in model_file:
                        # 提取第5列数据
                        columns = model_line.split()
                        if len(columns) > 4 and is_number(columns[4]):
                            with open(tK_file, 'a') as temp_file:
                                temp_file.write(columns[4] + '\n')

# 合并数据
tK_files = []
for t in temps:
    tK_file = f"{t}K"
    tK_files.append(tK_file)

paste_files(tK_files, "max_F_devi.dat")
