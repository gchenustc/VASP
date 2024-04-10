"""
Desc: 把二维材料卷成管子(x方向,向z轴上卷，卷完管子周期方向沿z轴，如果z轴堆叠多层，则只剩下一层), 第一个参数是结构文件名，第二到五个参数是预先对结构的扩胞倍数
Author: gchen
Time: 20240330
"""

# 通用的管子卷曲脚本

from ase.io import read,write
import sys
import numpy as np
from ase.build import make_supercell
import sys

def reverseAxis(atoms, n=2):
    """
		Reverse atoms along the corresponding axis (n=1:x, n=2:y, n=3:z), the reverse center is the cell center. 
	"""
    destAxisLen = np.linalg.norm(atoms.cell[n])
    for i in atoms:
        atomLen = np.dot(i.scaled_position, atoms.cell)[n]
        ratio = (destAxisLen-atomLen)/destAxisLen
        # print(i.position, atomLen, destAxisLen, atomLen/destAxisLen, ratio, i.position[2] * ratio)
        i.position[2] =  destAxisLen * ratio
    return atoms

def rollTube(name, sc, insideOut=True):
    """only support orthorhombic cell, and the tube is rolled along x to z direction."""
    atoms=read(name,format="vasp")
    atoms.center()
    
    if insideOut:
        reverseAxis(atoms, n=2)

    #* supercell
    P = np.array([[sc[0],0,0], [0,sc[1],0], [0,0,sc[2]]])
    atoms_sc = make_supercell(atoms, P)
    perimeter = atoms_sc.cell[0,0] # 沿着 x 轴卷曲后管子的周长
    R=perimeter/2/3.1415926 # 卷曲后的半径
    D_ang=2*np.pi/perimeter # 卷曲过程中应该旋转的弧度值
    center_z=atoms_sc.get_positions().mean(axis=0)[2]

    sum=0
    for i in atoms_sc:
        Ri=R-(i.z-center_z) # 原子在 x 轴上的偏移量
        theta=i.x*D_ang # 计算了每个原子在卷曲过程中应该旋转的角度 theta
        i.x=Ri*np.cos(theta)
        i.z=Ri*np.sin(theta)
        sum += 1
    print("Sum of Atoms =", sum)

    atoms_sc.set_cell([2*R+16, atoms_sc.cell[1,1], 2*R+16])
    atoms_sc.center()
    
    atoms_sc = atoms_sc[atoms_sc.numbers.argsort()]
    write("tube_"+name, atoms_sc)

name = sys.argv[1]
sc = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])]
rollTube(name, sc, insideOut=True)