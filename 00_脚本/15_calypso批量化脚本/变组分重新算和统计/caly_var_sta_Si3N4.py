"""
功能：获得基态结构，绘制convex_hull
首先要运行caly_var_verified.sh文件对calypso产生的文件重新优化。在results文件夹中运行该文件即可。
只支持单胞内含有两种元素的结构
此脚本是定制版，计算横坐标左侧是N和右侧是Si3N4的ratio
"""

from collections import Counter
import spglib as spg
import os
from ase.io import write,read
from pymatgen.io.vasp.outputs import Oszicar
from ase.phasediagram import PhaseDiagram
import numpy as np
import re
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt

### 基态结构分别为N和Si3N4 ###
### NxSiy = [N](x-y*(4/3))[N(4/3)Si](y) ### 中括号指的是formula，小括号是数量
# ----------- 编辑 -----------
symbols_list = ["N","Si"] # Si-->SiN(4/3)，这里Si不用换成SiN(4/3)
boundary_enth = {"N":-4.916493,"Si":-12.5717363}  # Si-->SiN(4/3),.这里Si不用换成SiN(4/3)，但能量给出的是Si3N4的能量
# ----------- end -----------

def get_sym(atoms, prec=0.1):
    """获得结构的对成性"""
    return spg.get_spacegroup(cell=atoms,symprec=prec)


def update_ehulls(strus):
    """更新ehull"""
    for atoms in strus:
        n_atoms = len(atoms.symbols)    # 原子数 [7,7,7,6,....]
        formula = atoms.info["formula"]  # 分子数，比如"N3P3"
        enthalpy_total = atoms.info["enthalpy"] * n_atoms
        # 获得atoms每种元素的数量 
        n_symbols = Counter(atoms.get_chemical_symbols()) # {"N":5, "P":4}

        # Si3N4(convex hull右侧基态)的比例
        # Si3N4 的比例: y/(x-(4/3)*y+y) = y/(x-(1/3)y)  # x指的是相对于纯N和纯Si基态来说N的比例，y指的是...纯Si的比例 
        y = n_symbols[symbols_list[1]]
        x = n_symbols[symbols_list[0]]
        ratio = y/(x-(1/3)*y)
        atoms.info["ratio"] = ratio

        # 基态是纯氮和纯Si: AxBy=xA+yB，ehull=(e(AxBy)-(x*e(A)+y*e(B)))/(x+y) 其中 minuend 就是 x*enthalpy(A)+y*enthalpy(B)
        # 基态是纯氮和SiN(4/3): ehull = e(NxSiy) - (y*e(SiN(4/3)) + (x-y*(4/3))*e(N))
        minuend = y*boundary_enth[symbols_list[1]] + (x-y*(4/3))*boundary_enth[symbols_list[0]]
        ehull = (enthalpy_total - minuend)/(x-(1/3)*y)
        atoms.info["ehull"] = ehull
        
def print_plot(strus, show=True):
    refs = []
    info = []
    for atoms in strus:
        n_atoms = len(atoms.symbols)
        formula = atoms.symbols.get_chemical_formula()
        ratio = atoms.info["ratio"]
        ehull = atoms.info["ehull"]

        spacegroup = atoms.info["spacegroup"]
        caly_index = atoms.info["caly_index"]
        index = atoms.info["index"]
        formula = atoms.info["formula"]
        enthalpy = atoms.info["enthalpy"]

        if 0<=ratio<=1:
            refs.append([ratio, ehull])
            info.append([ratio, formula, spacegroup, enthalpy, ehull, caly_index, index])
    
    # 打印信息
    print("id    ratio    formula    spacegroup      enthalpy     ehull              dir_path            index")
    info = sorted(info, key=lambda x: x[0])
    for idx,i in enumerate(info):
        print(f"{idx:<8d}{i[0]:<8.3f}{i[1]:<10s}{i[2]:<15s}{i[3]:<12f}{i[4]:<15f}{i[5]:<10s}\t{i[6]:<10s}")
    
    # 绘图
    plt.rcParams['lines.markersize'] = 25
    plt.rcParams['font.size'] = 40
    fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    figsize=(22, 16)
)

    points = np.array(refs)
    hull = ConvexHull(refs) #<scipy.spatial._qhull.ConvexHull object at 0x000002A5DC6D6820>
    ax.plot(points[:,0], points[:,1], 'o') # 先把所有的点画出来
    # 连线Hull上的点
    for simplex in hull.simplices:
        ratio = points[simplex, 0]
        e = points[simplex, 1]
            # simplex是二维坐标，[simplex, 0]指的是取points[simplex[0]]行和points[simplex[1]]]行合并后的第一列数据
        ax.plot(ratio, e, 'k-', lw=3)
    ax.set_xlabel(f"Si3N4 Ratio")
    ax.set_xlim(-0.01,1.01)
    plt.savefig("Convexhull.png", bbox_inches='tight', dpi=200)


def main(path, added_path, deep):
    strus_paths = []
    for root, dirs, files in os.walk(path):
        for _dir in dirs:
            if _dir in map(str, range(1,deep+1) ):
                contcar_path = os.path.join(os.path.join(root, _dir),'CONTCAR')
                oszicarcar_path = os.path.join(os.path.join(root, _dir),'OSZICAR')
                strus_paths.append((contcar_path,oszicarcar_path,root,_dir))

    for i in added_path:
        for root,dirs,files in os.walk(i):
            for _dir in dirs:
                contcar_path = os.path.join(os.path.join(root, _dir),'CONTCAR')
                oszicarcar_path = os.path.join(os.path.join(root, _dir),'OSZICAR')
                strus_paths.append((contcar_path,oszicarcar_path,root,_dir))

    strus = []
    for atoms in strus_paths:
        stru_ase = read(atoms[0], format="vasp")
        enthalpy_stru = Oszicar(filename=atoms[1]).final_energy
        
        stru_ase.info={}
        stru_ase.info["caly_index"]=atoms[2]
        stru_ase.info["index"]=atoms[3]
        stru_ase.info["formula"] = stru_ase.symbols.get_chemical_formula()
        stru_ase.info["enthalpy"]=float(enthalpy_stru)/len(stru_ase.get_chemical_symbols())
        stru_ase.info["spacegroup"]=get_sym(stru_ase)
        
        strus.append(stru_ase)

    return strus
    
path="./"  
added_path = ["supply_1","50GPa_1","50GPa_2","50GPa_3"]  # 额外添加的结构，需要在该文件下下再创建文件夹，里面是vasp计算的文件（需要CONTCAR和OSZICAR）
deep=5  # 验证caly_pso计算的深度，也就是每个比例下验证的数量
strus = main(path,added_path,deep) # 根据路径获得结构
update_ehulls(strus)
print_plot(strus)
