import os
from ase.io import read,write
from pymatgen.core import Structure
from ase.io import write,read
import numpy as np
import pymatgen.core as mg
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from pymatgen.analysis.structure_matcher import StructureMatcher



def _ase2pymatgen(ase_stru):
    """将ase的结构转换为pymatgen的结构"""
    if not os.path.isdir(".temp"):
        os.mkdir(".temp")
    random_name = f"{random.random()}.vasp"
    write(filename=f".temp/{random_name}", images=ase_stru, format="vasp")
    return Structure.from_file(filename=f".temp/{random_name}")

def _pymatgen2ase(pymatgen_stru):
    """将pymatgen的结构转换为ase的结构"""
    pymatgen_stru.to(filename=".temp.vasp", fmt="POSCAR")
    ase_stru = read(".temp.vasp")
    os.remove(".temp.vasp")
    return ase_stru
    
def symmetryCheck(atoms1, atoms2):
    # primitive_cell=True： 可以比较不同原子数的两个结构；scale：是否缩放到同一体积，高精度下关闭此项
    comp = StructureMatcher(ltol=1, stol=1, angle_tol=2, primitive_cell=True, scale=False)
    return comp.fit(atoms1, atoms2)