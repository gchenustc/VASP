import os
import numpy as np
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase import db
from ase.io import read, write
import pymatgen.core as mg
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
import random

def copyfile(path_1, path_2):
    """
    将路径path_1下的所有目录及文件复制到路径path_2下
    path_1: 待复制的目录或者文件路径
    path_2: 目标路径
    """
    if os.path.isdir(path_1):  # path_1是目录

        list_1 = os.listdir(path_1)
        if not list_1:  # 复制目录，仅仅是复制空目录
            os.mkdir(path_2)
        else:
            os.mkdir(path_2)  # 先复制最上层空目录
            for i in list_1:
                path_r = os.path.join(path_1, i)  # 下层目录或文件的绝对路径
                path_w = os.path.join(path_2, i)  # 目标目录或文件的绝对路径
                if os.path.isfile(path_r):  # 是文件则直接进行复制
                    with open(path_r, 'rb') as rstream:
                        container = rstream.read()
                        with open(path_w, 'wb') as wstream:
                            wstream.write(container)
                else:  # 是目录则调用本函数
                    copyfile(path_r, path_w)

    else:  # path_1是文件
        with open(path_1, 'rb') as rstream:  # 是文件直接复制文件
            container = rstream.read()
            file_name = os.path.basename(path_1)
            path_2 = os.path.join(path_2, file_name)
            with open(path_2, 'wb') as wstream:
                wstream.write(container)


def removedir(dir):
    dir = dir.replace('\\', '/')
    if(os.path.isdir(dir)):
        for p in os.listdir(dir):
            removedir(os.path.join(dir, p))
        if(os.path.exists(dir)):
            os.rmdir(dir)
    else:
        if(os.path.exists(dir)):
            os.remove(dir)
            
            
def duplicateCheck(db, atoms, check_range):
    """
    check_range: 传入db.select()中的参数，在其中进行重复检测
    """
    # 检测database是否为空，空则返回True
    if len(db) == 0:
        return False

    if check_range:
        atoms_to_be_check_list = list(
            map(lambda x: x.toatoms(), list(db.select(check_range))))
    else:
        atoms_to_be_check_list = list(
            map(lambda x: x.toatoms(), list(db.select())))

    return any(len(atoms.numbers) == len(dp_atoms_each.numbers) and np.all(atoms.numbers == dp_atoms_each.numbers) and np.all(np.abs(atoms.positions - dp_atoms_each.positions) < 0.001) for dp_atoms_each in atoms_to_be_check_list)


def ase2pymatgen(ase_stru):
    """将ase的结构转换为pymatgen的结构"""
    if not os.path.isdir(".temp"):
        os.mkdir(".temp")
    random_name = f"{random.random()}.vasp"
    write(filename=f".temp/{random_name}", images=ase_stru, format="vasp")
    return Structure.from_file(filename=f".temp/{random_name}")


def symmetryCheck(atoms_to_be_check_list, atoms):
    """
    atoms_to_be_check_list: 待检测的列表,pymatgen的格式
    atoms: 被检测的结构,pymatgen的格式
    如果atoms与atoms_to_be_check_list中的结构等价，返回True
    """
    if not atoms_to_be_check_list:
        return False
    comp = StructureMatcher()
    return any(comp.fit(ats, atoms) for ats in atoms_to_be_check_list)


def vaspMove(db_row, mode):
    if os.path.exists(f"./calc_record/id_{db_row.id}_{mode}"):
        removedir(f"./calc_record/id_{db_row.id}_{mode}")
    copyfile("workdir", f"./calc_record/id_{db_row.id}_{mode}")
    removedir("workdir")