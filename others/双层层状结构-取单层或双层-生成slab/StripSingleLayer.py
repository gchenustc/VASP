"""
作者：陈果
描述：取出两层层状结构的的单层，并加上20A的真空层（原子在底部，需要设置真空层所在轴）
"""
from pymatgen.core import  Structure as Sr
#import general_fun as gf 
from ase.io import write,read
from itertools import combinations
import os
import numpy as np
import random
import time
from sklearn.cluster import KMeans

def classfied(arr, class_=3, axis=1):
    # class_: 类别数
    # axis: 按照哪一个轴进行分类
    
    data = arr[:, axis]
    data = (data - np.mean(data)) / np.std(data)

    km = KMeans(n_clusters= class_)
    km.fit(data.reshape(-1,1))

    labels = km.labels_
    counts = np.bincount(labels) # 统计每个元素出现的次数
    min_indices = np.argsort(counts)[:2] # 找到出现次数较少的两个元素
    max_indice =  np.argsort(counts)[2]

    labels = np.where(np.isin(labels, min_indices), 100, labels) # 将出现次数较少的两个元素替换成指定的值
    labels = np.where(labels==max_indice, 1, labels) # 出现最多的元素替换成1
    labels = np.where(labels==100, 0, labels) # 出现次数较少的两个元素替换成0
    
    return labels

def write_vasp(stru, dis_between_upborder, axis, filename, labels, rmlabel):
    # 让一层之内的原子不分开（加了真空层后可能导致在周期性最上方和最下方的原子分离）
    scaled_pos = stru.get_scaled_positions()
    scaled_pos[:,axis] += dis_between_upborder + 0.01
    stru.set_scaled_positions(scaled_pos)
    stru.wrap(pbc=[True, True, True])  # 让周期性边界外的原子回到周期内

    pos = stru.get_positions() # 笛卡尔坐标
    
    rmindexs=[]
    for index, pos in enumerate(stru.positions):
        if labels[index]==rmlabel:
            rmindexs.append(index)
    for index in sorted(rmindexs,reverse=True):
        del stru[index]
    
    len_axis_vacuum = np.linalg.norm(stru.cell[axis])
    len_atoms_max = stru.positions.max()

    # ###########设置真空层为20A ###########
    add_vacuum = 20 - (len_axis_vacuum-len_atoms_max)      
    total_len_axis_vacuum = add_vacuum + len_axis_vacuum

    # stru.cell[1] 为真空层基矢，|x[a,b,c]| = total_len_axis_vacuum
    x = np.sqrt(total_len_axis_vacuum**2 / np.square(np.linalg.norm(stru.cell[axis])))
    stru.cell[axis] = x * stru.cell[axis]

    write(filename,stru, format="vasp")
    
if __name__ == "__main__":
    stru_path = "POSCAR"
    axis = 1 # 真空层所在轴
    stru_read = read(stru_path)

    # 获得层的分类
    labels = classfied(stru_read.get_scaled_positions(), class_=3, axis=axis)

    # 将为0的分类（有一半在最上方，另一半在最下方）向上移动到最下方的元素正好跨越到另一个周期，使得分类为0的元素不分开
    # 获得分类为0的元素的上半部分
    class0_up = []
    for i,pos in enumerate(stru_read.get_scaled_positions()):
        if labels[i]== 0 and pos[axis]>0.5:
            class0_up.append(pos.copy())
    dis_between_upborder0 = 1 - np.array(class0_up)[:,axis].min()


    # 将为1的分类向上移动到最下方的元素正好跨越到另一个周期
    # 获得分类为1的元素的上半部分
    class1_up = []
    for i,pos in enumerate(stru_read.get_scaled_positions()):
        if labels[i] == 1:
            class1_up.append(pos.copy())
    dis_between_upborder1 = 1 - np.array(class1_up)[:,axis ].min()

    stru_read_copy0 = stru_read.copy()
    write_vasp(stru_read_copy0, dis_between_upborder0, axis, "kind1.vasp", labels,1)

    stru_read_copy1 = stru_read.copy()
    write_vasp(stru_read_copy1, dis_between_upborder1, axis, "kind2.vasp", labels,0)