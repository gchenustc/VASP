import os
import shutil
import subprocess
"""
author: gchen
time: 20240523
desc: 将cp2k结构优化或者分子动力学产生的.pdb文件转换为.vasp文件，需要Multiwfn（执行命令是Multiwfn）
"""

with open(input("请输入要处理的文件名: "), 'r') as f:
    content = f.readlines()

# 指定要创建的文件夹路径
folder_path = "each_step_structures"

# 检查文件夹是否已经存在
if os.path.exists(folder_path):
    shutil.rmtree(folder_path)
    
if not os.path.exists(folder_path):
    try:
        # 创建文件夹
        os.makedirs(folder_path)
    except OSError as e:
        print(f"创建文件夹{folder_path}时出错: {e}")

split_lines = []
for line in content:
    if line.startswith("REMARK    Step"):
        split_lines.append(line)

start_index = int(input(f"请输入要从第几个 'REMARK    Step' 开始分割 (1-{len(split_lines)}): ")) - 1

for i in range(start_index, len(split_lines)):
    start = content.index(split_lines[i])  # 开始行的索引，从0开始
    
    if i == len(split_lines) - 1:
        end = len(content)
    else:
        end = content.index(split_lines[i+1]) - 1
    
    step_num = int(split_lines[i].split()[2][:-1])
    
    filename = folder_path+"/"+f"STEP_{step_num}.pdb"
    
    with open(filename, 'w') as out_file:
        out_file.write("TITLE pdb file\nAUTHOR gchen\n")
        out_file.writelines(content[start:end+1])
        
    print(f"已创建文件: {filename}")


with open("temp.sh", "w") as file:
    file.write("#!/bin/env bash\n\n")
    file.write(f"folder={folder_path}\n\n")
    file.write("cd $folder\n")
    file.write("for i in `ls`; do\n")
    file.write("    Multiwfn $i << EOF >/dev/null\n")
    file.write("    100\n")
    file.write("    2\n")
    file.write("    27\n")
    file.write("    ${i}.vasp\n")
    file.write("    0\n")
    file.write("    q\n")
    file.write("EOF\n")
    file.write("done\n")

subprocess.run(["bash", "temp.sh"], check=True)
