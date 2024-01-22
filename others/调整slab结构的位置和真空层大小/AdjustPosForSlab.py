"""
对于slab结构，调整原子位置，使之所有的原子都在底部（为了观察更方便），而且可以调整真空层厚度
"""
from pymatgen.core import  Structure as Sr
from ase.io import write,read
from itertools import combinations
import os
import numpy as np
import random
import time
from sklearn.cluster import KMeans

def classfied(arr, class_=2, axis=2):
    # class_: 类别数
    # axis: 按照哪一个轴进行分类
    
    data = arr[:, axis]
    data = (data - np.mean(data)) / np.std(data)

    km = KMeans(n_clusters= class_)
    km.fit(data.reshape(-1,1))

    labels = km.labels_
    
    return labels

def write_vasp(stru, dis_between_upborder=0, axis=2, filename="new.vasp", vacuum=20):
    # 让一层之内的原子不分开（加了真空层后可能导致在周期性最上方和最下方的原子分离）
    scaled_pos = stru.get_scaled_positions()
    scaled_pos[:,axis] += dis_between_upborder + 0.01
    stru.set_scaled_positions(scaled_pos)
    stru.wrap(pbc=[True, True, True])  # 让周期性边界外的原子回到周期内

    len_cell_axis = np.linalg.norm(stru.cell[axis])
    atoms_axis_max = stru.get_scaled_positions()[:,axis].max() * len_cell_axis
    cell_add = vacuum - (len_cell_axis - atoms_axis_max)
    total_len_cell_axis = len_cell_axis + cell_add

    # stru.cell[2] 为真空层基矢，|x[a,b,c]| = total_len_axis_vacuum
    x = np.sqrt(total_len_cell_axis**2 / np.square(np.linalg.norm(stru.cell[axis])))
    stru.cell[axis] = x * stru.cell[axis]

    write(filename,stru, format="vasp")
    
if __name__ == "__main__":
    stru_path = "new.vasp"
    axis = 2 # 真空层所在轴
    stru_read = read(stru_path)

    # 获得层的分类
    labels = classfied(stru_read.get_scaled_positions(), class_=2, axis=axis)

    # 判断在上层还是下层
    select1 = [i for i,j in enumerate(labels) if j ]
    stru_1 = stru_read.get_scaled_positions()[select1]
    select0 = [i for i,j in enumerate(1-labels) if j ]
    stru_0 = stru_read.get_scaled_positions()[select0]
    
    z_1_max = stru_1[:,axis].max()
    z_1_min = stru_1[:,axis].min()
    z_0_max = stru_0[:,axis].max()
    z_0_min = stru_0[:,axis].min()
 
    if z_1_min > z_0_max and (z_1_min - z_0_max)*np.linalg.norm(stru_read.cell[axis]) > 3:
        up_arr = stru_1
        flag = 0
    elif z_0_min > z_1_max and (z_0_min - z_1_max)*np.linalg.norm(stru_read.cell[axis]) > 3:
        up_arr = stru_0
        flag = 0
    else:
        up_arr = stru_read.get_scaled_positions()
        flag = 1
#    print(up_arr)
    
    # 如果分不开两类，则向下移动
    if flag:
        dis_between_upborder = -np.array(up_arr)[:,axis].min()
    else:
        dis_between_upborder = 1 - np.array(up_arr)[:,axis ].min()
    write_vasp(stru_read, dis_between_upborder, axis, filename="new1.vasp", vacuum=10)  # 调整真空层大小