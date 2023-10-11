from collections import Counter
import spglib as spg
import os
from ase.io import write,read
from pymatgen.io.vasp.outputs import Oszicar
from ase.phasediagram import PhaseDiagram
import numpy as np
import re


# ----------- 编辑 -----------
symbols_list = ["N","Si"]
# ----------- end -----------

def get_sym(atoms, prec=0.1):
    """获得结构的对成性"""
    ret = spg.get_spacegroup(cell=atoms,symprec=prec)
    return "None" if ret==None else ret 


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

        
def print_plot(strus, show=True):
    refs = []  # eg: [("N3",ehull1),("N3P3",ehull2)]
    info = []
    for atoms in strus:
        formula = atoms.info["formula"]  # 分子数，比如"N3P3"
        n_atoms = len(atoms.symbols)    # 原子数 [7,7,7,6,....]
        enthalpy_total = atoms.info["enthalpy"] * n_atoms
        spacegroup = atoms.info["spacegroup"]
        caly_index = atoms.info["caly_index"]
        index = atoms.info["index"]
        enthalpy = atoms.info["enthalpy"]

        info.append((formula, spacegroup, enthalpy, caly_index, index))
        info = sorted(info, key=lambda x: x[2])
        #print(formula, spacegroup, enthalpy, caly_index, index)

    
    # 打印信息
    for idx,i in enumerate(info):
        print(f"{idx:<8d}{i[0]:<10s}{i[1]:<15s}{i[2]:<12f}{i[3]:<10s}\t{i[4]:<10s}")
 


strus_paths = []
for root, dirs, files in os.walk("./"):
    for _dir in dirs:
        if _dir in [str(i) for i in range(1,51)]:
            contcar_path = os.path.join(os.path.join(root, _dir),'CONTCAR')
            oszicarcar_path = os.path.join(os.path.join(root, _dir),'OSZICAR')
            strus_paths.append((contcar_path,oszicarcar_path,root,_dir))

strus = []
for atoms in strus_paths:
    try:
        stru_ase = read(atoms[0], format="vasp")
    except (ValueError,IndexError):
        continue
    enthalpy_stru = Oszicar(filename=atoms[1]).final_energy
    
    stru_ase.info={}
    stru_ase.info["caly_index"]=atoms[2]
    stru_ase.info["index"]=atoms[3]
    stru_ase.info["formula"] = stru_ase.symbols.get_chemical_formula()
    stru_ase.info["enthalpy"]=float(enthalpy_stru)/len(stru_ase.get_chemical_symbols())
    stru_ase.info["spacegroup"]=get_sym(stru_ase, 0.1)
    #print(stru_ase.info["spacegroup"])
    
    strus.append(stru_ase)
    
strus = remove_repeated_strus(strus)
print_plot(strus)
