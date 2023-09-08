from ase.io import write,read
import numpy as np
import pymatgen.core as mg
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
import os
import argparse
"""
description: 检测两个结构是否为等价，支持多种文件格式
usage:
python *.py -c1 POSCAR1 -c2 POSCAR2  # POSCAR1和POSCAR2分别分两个结构文件的路径
"""

parser = argparse.ArgumentParser(description="None")
parser.add_argument('-c1','--POSCAR1', default=["POSCAR1"], help="Specify the POSCAR1 file")
parser.add_argument('-c2','--POSCAR2', default=["POSCAR2"], help="Specify the POSCAR2 file")
prm = parser.parse_args()

def _ase2pymatgen(ase_stru):
    """将ase的结构转换为pymatgen的结构"""
    write(filename=".vasp", images=ase_stru, format="vasp")
    pymatgen_stru = Structure.from_file(filename=".vasp")
    os.remove(".vasp")
    return pymatgen_stru

def symmetryCheck(atoms1, atoms2):
    
    atoms1 = _ase2pymatgen(read(atoms1))
    atoms2 = _ase2pymatgen(read(atoms2))

    comp = StructureMatcher()
    return comp.fit(atoms1, atoms2)
    
print("是否为相同结构：",symmetryCheck(prm.POSCAR1, prm.POSCAR2))