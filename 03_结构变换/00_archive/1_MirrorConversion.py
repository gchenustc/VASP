"""
作用：
main('POSCAR','z') 指的是分别沿着对z轴方向进行中心对称
注意：
只支持 vasp 的POSCAR 格式
且POSCAR是笛卡尔坐标
"""

import numpy as np
import os
from math import sqrt,pow
import shutil

def retPoscar(poscar, const, pos):
    """
    输出指定的 POSCAR
    """
    atomsNum = pos.shape[0]
    
    with open(poscar,'r') as f:
        head = f.readlines()
        head_copy = head[:] # bak

        head = head[:8]
        head[2] = '    '.join(str(const[0].tolist()).strip('[]').split(',')) + '\n'
        head[3] = '    '.join(str(const[1].tolist()).strip('[]').split(',')) + '\n'
        head[4] = '    '.join(str(const[2].tolist()).strip('[]').split(',')) + '\n'
        head[6] = str(atomsNum) + '\n'
        with open('b.vasp.new','w') as fw:
            fw.writelines(head)
            np.savetxt(fw,pos,fmt='%.8f')


def mirrorConversion(pos, const, axis):
    
    value_ranges = ['x','X','y','Y','z','Z']
    axis_index_list = [0,0,1,1,2,2]
    
    assert axis in value_ranges
    
    corres_dict = dict((zip(value_ranges,axis_index_list)))
    
    axis_index = corres_dict[axis]
    const_value = const[axis_index][axis_index]
    pos[:,axis_index] = const_value - pos[:,axis_index]
    

def main(poscar,axis='x'):

    pos_ori = np.loadtxt(poscar,dtype=np.float64,skiprows=8)
    const = np.loadtxt(poscar,dtype=np.float64,skiprows=2,max_rows=3)

    pos = pos_ori.copy()
    mirrorConversion(pos, const, axis)

    retPoscar(poscar, const, pos)
    
main('POSCAR','z')