import numpy as np
from pymatgen.io.vasp.outputs import Xdatcar
import os,shutil
import time

"""
Author: gchen
Time: 2024.9.13
Desc: This script is used to eliminate centroid offset of XDATCAR file.
"""
INPUT_FILE_NAME = "XDATCAR"
OUTPUT_FILE_NAME = "XDATCAR_o"

def rm_posiline_element(lines):
    """remove the element name in the posiline. eg. 0.000000 0.000000 0.000000 N, remove N"""
    new_lines = []  
    for line in lines:  
        if len(line.split()) ==4  and line.split()[0].replace('.', '', 1).replace('-', '').isdigit():  # 判断是否为原子位置行  
            # 去掉元素名称  
            new_lines.append(' '.join(line.split()[:-1]) + '\n')  
        else:  
            new_lines.append(line) 
    return new_lines


def change_line(file_path, n_line, new_line_content): 
    
    """chage the n_line of file_path to new_line_content"""
    
    # 读取文件内容  
    with open(file_path, 'r', encoding='utf-8') as file:  
        lines = file.readlines()  # 读取所有行到一个列表中  

    if len(lines) >= n_line:  # 确保文件至少有n_line行  
        lines[n_line-1] = new_line_content  # 替换内容  
    
    lines = rm_posiline_element(lines)
    
    return lines


def cont_file(source_file, target_file):

    """copy the content of source_file to target_file"""
    
    with open(source_file, 'r') as src:  # 以读模式打开源文件  
        content = src.read()  # 读取源文件内容  

    with open(target_file, 'a') as target:  # 以追加模式打开目标文件  
        target.write(content)  # 将内容写入目标文件 


def xdatcar_operate(xdatcar):
    
    """operate the xdatcar"""
    
    ref_pos = xdatcar.structures[0].frac_coords # reference position: the 0 frame, first atom

    for i in range(len(xdatcar.structures)):
    # if i > 2: break
    # find the bias of the reference position
        bias = (ref_pos - xdatcar.structures[i].frac_coords).mean(axis=0)

    # frac_bias = xdatcar.structures[i].lattice.get_fractional_coords(bias)

        for j in range(len(xdatcar.structures[i])):
            element = xdatcar.structures[i][j].specie.value
            fixed_coord = xdatcar.structures[i].frac_coords[j] + bias
            xdatcar.structures[i][j] = element, fixed_coord

if __name__ == '__main__':
    
    start_time = time.perf_counter()
    
    xdatcar = Xdatcar(INPUT_FILE_NAME)

    xdatcar_operate(xdatcar)
    print("*** Reading XDATCAR file successfully.")

    tmp_folder_name = ".tmp.XDATCAR"

    os.mkdir(tmp_folder_name) 

    tmp_file_name = []

    # write the structures to the tmp folder
    for i in range(len(xdatcar.structures)):
        # if i > 2: break
        
        file_name = os.path.join(tmp_folder_name, f".{i+1}.vasp")
        tmp_file_name.append(file_name)
        xdatcar.structures[i].to(fmt='poscar', filename=file_name)
    print("*** Operating XDATCAR file successfully.")

    # create a new file for new XDATCAR.
    with open(OUTPUT_FILE_NAME, 'w') as file:  
        pass
    
    # change the "Direct configuration=     n" in the XDATCAR file of line 8 for each step.
    for i in range(len(tmp_file_name)):
        lines = change_line(tmp_file_name[i], 8, f"Direct configuration=     {i+1}\n")
            
        with open(tmp_file_name[i], 'w', encoding='utf-8') as file:  
            file.writelines(lines)
            
        cont_file(tmp_file_name[i], OUTPUT_FILE_NAME)
    print("*** Writing XDATCAR file successfully.")

    # remove the tmp folder
    try:
        shutil.rmtree(tmp_folder_name) 
    except PermissionError as e:
        print(f"PermissionError: {e}")
    print("*** Removing the tmp folder successfully.")
        
    time_cost = time.perf_counter() - start_time
    print(
        f'*** Time used: {time.strftime("%H:%M:%S", time.gmtime(time_cost))}'
    )

