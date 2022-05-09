# -*- coding:utf-8 -*-
"""
 author：gchen 
 time:2021.11.04.20:24:01
"""
import numpy as np
import pandas as pd
import os
import math

def get_path():
    """获得文件路径"""
    file_name = input("请输入文件名:")
    path = os.getcwd() + os.sep + file_name
    return path

def get_head(path):
    """获取文件头"""
    file_head = []
    with open(path,"r") as f:
        count=8
        while count>0:
            count -= 1
            file_head.append(f.readline())
    return file_head

def num_adjust(num):
    """对每一个原子坐标进行相同长度的切割（即保证对其）,输出每一个坐标信息，以列表形式保存"""
    num_tran = str(num).ljust(25,"0")[0:10]
    date_out.append(num_tran)

def num_adjust_temp(num):
    num_tran = str(num).ljust(25,"0")[0:10]
    date_out_temp.append(num_tran)

def to_skip_atoms(skip_atoms):
    date_surface_atoms.extend(date_out_temp[len(date_out_temp)-skip_atoms*3:len(date_out_temp)])

def sort_z(path,atoms):
    """对z轴进行排序,返回完整坐标数据和原子数"""
    data = pd.read_table(path,skiprows=8,delimiter="\s+\s+",header=None,engine="python",\
                         names=["x","y","z"],dtype=np.float64)

    (data.T).applymap(num_adjust_temp)

    data_copy = data.copy()
    for i in range(atoms):
        data_copy.drop(list(data_copy.index)[-1],inplace=True)
    data_copy.sort_values("z",inplace=True,ascending=False)
    ###输入坐标信息###
    (data_copy.T).applymap(num_adjust)

    data.sort_values("z", inplace=True, ascending=False)

    #以txt形式输出(临时文件)
    path_ = os.path.split(path)[0] + os.sep + "temp.txt"
    data.to_csv(path_,header=False,index=False,sep=" ")
    #读取坐标文件(从临时文件中读取)
    with open(path_,"r") as f:
        data = []  #data是原子坐标数据，初始化，再一行一行读取
        atoms = 0  #vasp文件原子数目，每读取一行，atoms加1
        while True:
            line = f.readline()
            if line == "": break
            atoms += 1
            line_new = line.replace(" ","    ")
            data.append(line_new)
    #删除临时文件
    os.remove(path_)
    return data,atoms

def set_layers(atoms,skip_atoms):
    """对原子进行分层,返回(层数，每层的原子数)元组"""
    layers_rec = dict()
    for i in range(2,(atoms-skip_atoms)//4):
        num = (atoms-skip_atoms)/i
        if num%1 == 0:
            layers_rec[i] = int(num)
    print("一共有%d个原子，参与固定的原子有%d个，推荐分为%r层，每层%r个原子" % (atoms,atoms-skip_atoms,list(layers_rec.keys()),list(layers_rec.values())))

    while True:
        layers = float(input("请输入需要分的层数:"))
        print("您分了%.2f层，每层%.4f个原子。" % (layers, (atoms-skip_atoms)/layers))
        if input("确定这样分吗？(y)").lower() == "y": break

    return layers, (atoms-skip_atoms)/layers

def layers_fixed():
    """确定固定原子的层数"""
    return int(input("请输入需要固定的层数:"))

def fix_atoms(data,atoms_begin_fix,atoms_end_fix):
    global POSCAR

    atom_counts = 1
    for i in range(1,len(data)+1):
        POSCAR = POSCAR + data[i-1] + "    "
        if not i%3:
            if atoms_begin_fix<= atom_counts <= atoms_end_fix:
                POSCAR += "    F   F   F"
            else:
                POSCAR += "    T   T   T"
            atom_counts += 1
            POSCAR += "\n"

def main():
    print("-----推荐使用笛卡尔坐标-----")
    path = get_path()         #POSCAR文件路径
    file_head = get_head(path)  #POSCAR的前八行

    global POSCAR
    skip_switch = input("分层要跳过表面原子吗(y or n):").lower()
    if skip_switch == "y":
        skip_atoms = int(input("请输入要跳过的表面原子数:"))
    else:
        skip_atoms = 0
    data,atoms = sort_z(path,skip_atoms)   #(1)原子坐标信息,以列表形式储存(2)POSCAR内的原子数
    if skip_atoms:
        to_skip_atoms(skip_atoms)

    layers,layer_atoms = set_layers(atoms,skip_atoms)   #层数以及每层的原子数
    num_layer_fix = layers_fixed()     #固定的层数
    num_atoms_fix = math.floor(num_layer_fix*layer_atoms)  #固定的原子数

    #如果分的层数是小数的话，则做修正
    if layers%1:
        num_atoms_fix +=1
    """
    coor_begin_fix = (atoms - num_atoms_fix)//2 + 1  #确定从第n个坐标开始固定
    coor_end_fix = coor_begin_fix + num_atoms_fix -1 #确定从第n个坐标结束固定
    """
    #从中间固定
    atoms_begin_fix_middle = ((layers - num_layer_fix)/2)*layer_atoms + 1  #确定从第n个原子开始固定
    atoms_end_fix_middle = atoms_begin_fix_middle + num_atoms_fix -1 #确定第n个原子结束固定
    #从底部固定
    atoms_begin_fix_bottom = atoms-skip_atoms - num_atoms_fix + 1  #确定从第n个原子开始固定
    atoms_end_fix_bottom =  atoms-skip_atoms   #确定第n个原子结束固定

    POSCAR=""
    for i in range(len(file_head)-1):
        POSCAR += file_head[i]
    POSCAR += "Selective dynamics\n"
    POSCAR += file_head[-1]

    fixed_way = input("请输入固定的方式(a：底部固定，b：中间固定):")
    """
    方法2：
    if fixed_way.lower() == "a":
        begin = atoms_begin_fix_bottom
        end = atoms_end_fix_bottom
    elif fixed_way.lower() == "b":
        begin = atoms_begin_fix_middle
        end = atoms_end_fix_middle
    else:
        raise TypeError

    
    for i in range(1,len(data)+1):
        line_data = data[i-1].replace(" ","    ")
        line_data = line_data[:-2].ljust(45," ") + "\n"
        if begin <= i <=end:
            line_data = line_data.replace("\n","    F  F  F\n")
        POSCAR += line_data    
    """
    if fixed_way.lower() == "b":
        fix_atoms(date_out, atoms_begin_fix_middle, atoms_end_fix_middle)
    elif fixed_way.lower() == "a":
        fix_atoms(date_out, atoms_begin_fix_bottom, atoms_end_fix_bottom)
    else: raise TypeError
    if skip_atoms:
        fix_atoms(date_surface_atoms,0,0)

    #输出POSCAR
    with open(os.path.split(path)[0] +os.sep +"POSCAR","w") as f:
        f.write(POSCAR)

if __name__ == "__main__":
    date_out=[]
    date_out_temp = []
    date_surface_atoms = []
    main()




