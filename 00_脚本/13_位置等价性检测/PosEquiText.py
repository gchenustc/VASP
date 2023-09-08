"""
author: gchen
date: 2023.9.8
desc: 原子位点等价性检测，命令行中直接运行该文件
需要修改 POSCAR文件名和指定检测的某种元素类别
"""
from ase.io import write,read
import numpy as np
import pymatgen.core as mg
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
import os
import itertools as it

# --- 修改 ---
poscar = read("POSCAR",format="vasp")  ## POSCAR文件名
select_kind = "all"  ## 需要对比的原子类别，"all" 检测所有
# END

def _ase2pymatgen(ase_stru):
    """将ase的结构转换为pymatgen的结构"""
    write(filename=".vasp", images=ase_stru, format="vasp")
    pymatgen_stru = Structure.from_file(filename=".vasp")
    os.remove(".vasp")
    return pymatgen_stru

def symmetryCheck(atoms1, atoms2):
    atoms1 = _ase2pymatgen(atoms1)
    atoms2 = _ase2pymatgen(atoms2)
    
    # primitive_cell=True： 可以比较不同原子数的两个结构；scale：是否缩放到同一体积，高精度下关闭此项
    comp = StructureMatcher(ltol=0.2, stol=0.2, angle_tol=2, primitive_cell=True, scale=False)
    return comp.fit(atoms1, atoms2)

select_idx = []
for idx,atom in enumerate(poscar):
    if select_kind.strip() == "all":
        select_idx.append(idx)
    else:
        if atom.symbol == select_kind:
            select_idx.append(idx)
        

combi = it.combinations(select_idx,2)

for idx1,idx2 in combi:
    poscar1 = poscar.copy()
    poscar2 = poscar.copy()
    poscar1[idx1].symbol = "Cl"
    poscar2[idx2].symbol = "Cl"
    print(f"no.{idx1+1} ({' '.join(map(str,poscar1[idx1].position))})\n-\n\
no.{idx2+1} ({' '.join(map(str,poscar2[idx2].position))}); \n-\n是否等价：{symmetryCheck(poscar1,poscar2)}")
    print("-"*10)


