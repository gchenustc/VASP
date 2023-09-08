"""
author: gchen
date: 2023.9.8
desc: 指定traj文件路径，和搜索结构的种类（最多两种），统计这些结构并生成convexhull，保存基态线上的结构，
注意在traj文件中需要包含boundary的结构，用于计算ehull
"""
import pickle
from collections import Counter
from ase.io import Trajectory
import spglib as spg
import os
from ase.io import write
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from ase.phasediagram import PhaseDiagram
import pickle
import numpy as np


### *       输入参数         !###
#! 指定要遍历的文件夹路径
folder_path = './'
#! 识别对称性的精度
prec = 0.2
#! 包含的原子种类列表(有且只有两种元素)
symbols_list = ["N", "P"]
#! 打开的储存的结构的文件类型
suffix = "traj"
### *       输入参数         !###


class GetDestFiletypePaths(object):
    """获得指定文件类型的所有路径"""

    def __init__(self, folder_path):
        self.folder_path = folder_path  # 指定文件夹路径
        self.file_paths = []  # 只要是该路径下的文件路径都出存在着
        self.dest_file_paths = []  # 只包含指定后缀的文件路径

    def get_file_paths(self):
        """输入文件夹，获得该文件夹以及子文件夹中的所有文件"""
        # * 获得指定目录下的所有文件
        # * 使用walk函数遍历文件夹中的文件
        # * root: 目录, dirs: 目录下的文件夹, file: 目录下的文件名
        for root, dirs, files in os.walk(self.folder_path):
            for file in files:
                file_path = os.path.join(root, file)
                self.file_paths.append(file_path)

    def clear_file_paths(self):
        """清除所有文件的路径"""
        self.file_paths = []
        self.dest_file_paths = []  # 只包含指定后缀的文件路径

    def judge_file_type(self, path, suffix):
        """
        确定目标路径是否为".suffix"结构

        Parameters
        ----------
        path - str
            文件路径
        suffix - str
            需要判断的后缀类型

        Return
        ----------
        True or False
        """
        suffix_file = path.split(".")[-1]
        return suffix_file == suffix

    def get_dest_filetype_paths(self, suffix):
        """输入文件后缀，获得该文件后缀的所有文件列表"""
        self.get_file_paths()  # 获得所有文件的路径
        self.dest_file_paths.extend([
            path for path in self.file_paths if self.judge_file_type(path, suffix)])


class GetEhull(object):
    """获得ehull"""

    def __init__(self, prec, symbols_list):
        self.prec = prec
        self.symbols_list = symbols_list
        # * boundary结构的ehull值，以及对应的atoms结构，后面会更新
        # * eg: boundary_ehth = {"N":1000; "P":1000}
        self.boundary_enth = {i: 1000 for i in self.symbols_list}
        # * eg: boundary_enth_strus = {"N": <type 'atoms1'>, "P": <type 'atoms2'>}
        self.boundary_enth_strus = {i: None for i in self.symbols_list}
        # ! 去除了P1和重复结构的列表 - 最终要的结果
        self.list_prop_strus = []

    def remove_repeated_strus(self):
        """去除重复的结构
        Parameters
        ----------
        strus_list - list of ase.atoms
            存放结构的列表

        Return
        ----------
        返回去重的结构列表
        """
        ret = []
        for stru in self.list_prop_strus:
            if stru not in ret:
                ret.append(stru)
        return ret

    def updata_boundary(self, atoms):
        """更新边界（单一组分）的结构和焓值信息，该信息用来计算ehull"""
        symbols = atoms.get_chemical_symbols()  # ["P","N","N",...]
        if len(set(symbols)) == 1:
            symbol = symbols[0]
            if atoms.info["enthalpy"] < self.boundary_enth[symbol]:
                self.boundary_enth[symbol] = atoms.info["enthalpy"]
                self.boundary_enth_strus[symbol] = atoms

    def collect_dest_atoms(self, file_paths, openway=Trajectory, *args):
        """
        收集指定类型的结构

        Parameters
        ----------
        trajs_path - list of str
            文件的路径

        openway - function
            打开file的方式
        """

        for file in file_paths:
            atoms_list = openway(file, *args)  # 打开file
            for atoms in atoms_list:
                # ! 判断该结构的类型是否为P1
                if spg.get_symmetry_dataset(atoms, symprec=self.prec)['number'] != 1:
                    self.list_prop_strus.append(atoms)

                    # * 获得当前非P1结构的组分，如果原子类型不超过1，则可能是boundary结构，则更新boundary结构的信息
                    self.updata_boundary(atoms)
        self.remove_repeated_strus()

    def update_ehulls(self):
        """更新ehull"""
        for atoms in self.list_prop_strus:
            n_atoms = len(atoms.symbols)    # 原子数
            formula = atoms.symbols.get_chemical_formula()  # 分子数，比如"N3P3"
            enthalpy_total = atoms.info["enthalpy"] * n_atoms
            # 获得atoms每种元素的数量 {"N":5, "P":4}
            n_symbols = Counter(atoms.get_chemical_symbols())

            # AxBy=xA+yB，ehull=(enthalpy(AxBy)-(x*enthalpy(A)+y*enthalpy(B)))/(x+y)
            # 其中 minuend 就是 x*enthalpy(A)+y*enthalpy(B)
            minuend = n_symbols[self.symbols_list[0]] * self.boundary_enth[self.symbols_list[0]
                                                                           ] + n_symbols[self.symbols_list[1]] * self.boundary_enth[self.symbols_list[1]]
            ehull = (enthalpy_total - minuend)/n_atoms
            atoms.info["ehull"] = ehull

    def plot_convex_hull(self, show=True):
        refs = []  # eg: [("N3",ehull1),("N3P3",ehull2)]
        for atoms in self.list_prop_strus:
            formula = atoms.symbols.get_chemical_formula()
            ehull = atoms.info["ehull"]
            refs.append((formula, ehull))
        self.pd = PhaseDiagram(refs)
        self.pd.plot(show=show)

    def save_ehull_strus_to_vasp(self):
        """保存在ehull线上的结构"""
        index = np.where(
            self.pd.hull == True)  # 判断结构是否在ehull线上 eg: [[True, False, ....]]
        for i in index[0]:
            write(f"{i}.vasp", self.list_prop_strus[i], format="vasp")
            print("info:", i, self.list_prop_strus[i].info["enthalpy"], spg.get_spacegroup(
                self.list_prop_strus[i], symprec=self.prec))

    def save_strus(self):
        with open("out.pkl", "wb") as f:
            pickle.dump(self.list_prop_atoms, f)


if __name__ == "__main__":
    cls_path = GetDestFiletypePaths(folder_path)
    cls_ehull = GetEhull(prec, symbols_list)

    # * 获得指定的file_type的所有文件路径
    cls_path.get_dest_filetype_paths(suffix)

    # * 收集指定的结构(非P1的结构)
    # 指定文件路径，打开该文件的函数，该函数的其他参数-"r"
    cls_ehull.collect_dest_atoms(cls_path.dest_file_paths, Trajectory, "r")

    # * 更新ehull参数
    cls_ehull.update_ehulls()

    # * 绘制convex hull
    cls_ehull.plot_convex_hull()

    # * 保存在ehull上的结构，vasp类型
    cls_ehull.save_ehull_strus_to_vasp()
