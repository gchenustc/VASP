"""
获得基态结构，绘制convex_hull
首先要运行caly_var_recalc.sh文件对calypso产生的文件重新优化。此文件和caly_var_recalc.sh在一个文件夹内
只支持单胞内含有两种元素的结构
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

# ----------- 编辑 -----------
symbols_list = ["N","Si"]
boundary_enth = {"N":-4.916493,"Si":-1.209191}
# ----------- end -----------

def get_sym(atoms, prec=0.1):
    """获得结构的对成性"""
    return spg.get_spacegroup(cell=atoms,symprec=prec)


def remove_repeated_strus(strus_list):
    """去除重复的结构
    Parameters
    ----------
    strus_list - list of ase.atoms
        存放结构的列表

    """
    ret = []
    for stru in strus_list:
        if stru not in ret:
            ret.append(stru)
    return ret


def update_ehulls(strus):
    """更新ehull"""
    for atoms in strus:
        n_atoms = len(atoms.symbols)    # 原子数 [7,7,7,6,....]
        formula = atoms.info["formula"]  # 分子数，比如"N3P3"
        enthalpy_total = atoms.info["enthalpy"] * n_atoms
        # 获得atoms每种元素的数量 
        n_symbols = Counter(atoms.get_chemical_symbols()) # {"N":5, "P":4}
        
        # 元素含量比例
        ratio = n_symbols[symbols_list[1]] / (n_symbols[symbols_list[1]] + n_symbols[symbols_list[0]])
        atoms.info["ratio"] = ratio

        # AxBy=xA+yB，ehull=(enthalpy(AxBy)-(x*enthalpy(A)+y*enthalpy(B)))/(x+y)
        # 其中 minuend 就是 x*enthalpy(A)+y*enthalpy(B)
        minuend = n_symbols[symbols_list[0]] * boundary_enth[symbols_list[0]
                                                                       ] + n_symbols[symbols_list[1]] * boundary_enth[symbols_list[1]]
        ehull = enthalpy_total - minuend
        atoms.info["ehull"] = ehull

        
def print_plot(strus, show=True):
    refs_pd = []  # eg: [("N3",ehull1),("N3P3",ehull2)]
    refs = []
    info = []
    for atoms in strus:
        n_atoms = len(atoms.symbols)
        formula = atoms.symbols.get_chemical_formula()
        ehull = atoms.info["ehull"]
        ratio = atoms.info["ratio"]
        refs.append([ratio, ehull/n_atoms])
        refs_pd.append([formula,ehull])
        
        spacegroup = atoms.info["spacegroup"]
        caly_index = atoms.info["caly_index"]
        index = atoms.info["index"]
        formula = atoms.info["formula"]
        enthalpy = atoms.info["enthalpy"]
        

        info.append([ratio, formula, spacegroup, enthalpy, ehull/n_atoms, caly_index, index, atoms])
    
    # 打印信息
    pd = PhaseDiagram(refs_pd,verbose=False)
    ehull_index = np.where(pd.hull == True)[0]   
    for i,j in enumerate(info):
        if i in ehull_index:
            j.append(True)
        else:
            j.append(False)
        
    print("id    ratio    formula    spacegroup      enthalpy     ehull              dir_path            index        if_in_hull")
    info = sorted(info, key=lambda x: x[0])
    for idx,i in enumerate(info):
        print(f"{idx:<8d}{i[0]:<8.3f}{i[1]:<10s}{i[2]:<15s}{i[3]:<12f}{i[4]:<15f}{i[5]:<10s}\t{i[6]:<10s}{i[8]}")
    
    # 绘图
    plt.rcParams['lines.markersize'] = 20
    plt.rcParams['font.size'] = 35
    fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    figsize=(22, 16)
)
    
    points = np.array(refs)
    hull = ConvexHull(refs)
    ax.plot(points[:,0], points[:,1], 'o')
    for simplex in hull.simplices:
        ax.plot(points[simplex, 0], points[simplex, 1], 'k-')
    #print(hull.good)
    ax.set_xlabel(f"{symbols_list[1]} ratio")
    plt.savefig("convexhull.png", bbox_inches='tight', dpi=200)    
    
    # 保存ehull结构
    for idx,i in enumerate(info):
        if i[8]:
            write(f"{idx}.vasp", i[7], format="vasp")
        

def main(path, added_path, deep):
    # 统计额外添加的文件
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
added_path = ["50GPa_1"]  # 额外添加的结构，需要在该文件下下再创建文件夹，里面是vasp计算的文件（需要CONTCAR和OSZICAR）
deep=5  # 验证caly_pso计算的深度，也就是每个比例下验证的数量
strus = main(path,added_path,deep) # 根据路径获得结构
#strus = remove_repeated_strus(strus)
update_ehulls(strus)
print_plot(strus)

