"""
作用：main('POSCAR',x,y,z) 指的是对POSCAR分别沿着a,b,c轴移动x,y,z埃

只支持 vasp 的 POSCAR 格式
且 POSCAR是笛卡尔坐标
"""

import numpy as np
import os
from math import sqrt,pow
import shutil

def retPoscar(poscar, const, arr):
    """
    输出指定的 POSCAR
    """
    atomsNum = arr.shape[0]
    
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
            np.savetxt(fw,arr,fmt='%.8f')

def main(poscar, upmove=0, forwardmove=0, rightmove=0):

    arr_ori = np.loadtxt(poscar,dtype=np.float64,skiprows=8)
    const = np.loadtxt(poscar,dtype=np.float64,skiprows=2,max_rows=3)

    arr = arr_ori.copy()
    arr[:,0] += upmove
    arr[:,1] += forwardmove
    arr[:,2] += rightmove

    retPoscar(poscar, const, arr)

if __name__ == "__main__":
    main('POSCAR',0,0,0)