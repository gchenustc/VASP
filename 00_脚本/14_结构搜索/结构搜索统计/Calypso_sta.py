"""
author: gchen
date: 2023.9.8
desc: 指定Calypso储存结构的文件路径，搜索结构的种类（最多两种），以及boundary的能量，指定的内容在脚本的末尾更改。统计这些结构并生成convexhull，保存基态线上的结构。
注意Calypso储存结构的文件不能包含boundary的结构，boundary能量需要手动输入
"""
from collections import Counter
import spglib as spg
import os
from ase.io import write,read
from ase.phasediagram import PhaseDiagram
import numpy as np
import re

class File2class(object):
    """将指定文件夹中的calypso结构文本转换为list[ase.atoms]"""
    def __init__(self, folder_path="./", element=["N","H"], prec=0.1):
        self.folder_path = folder_path
        #
        self.element=element
        self.prec=prec
        
    def calc_path(self):
        """获得指定文件夹中所有文件的路径"""
        self.file_paths=[ ]
        for root, dirs, files in os.walk(self.folder_path):
            for file in files:
                file_path = os.path.join(root, file)
                self.file_paths.append(file_path)
    
    def calc_dest_path(self, re_pattern=".*pso_opt_[\d]+$"):
        """获得指定文件夹中指定文件（calypso储存结构的文件）的路径"""
        self.dest_file_paths=[]
        for file in self.file_paths:  
            match = re.match(re_pattern, file)
            if match != None:
                self.dest_file_paths.append(match.group())     
    
    def calc_atoms_dict(self):
        """分割calypso结构文件，提取POSCAR文本"""
        self.atoms_dict = {} # 文件名：list[poscar结构文本]
        self.enthalpy_dict = {}  #文件名：list[poscar焓值]
        
        for file in self.dest_file_paths:
            
            key_name = file.split("/")[-1] # 文件名，去掉文件名中的./；
            self.atoms_dict[key_name] = []
            self.enthalpy_dict[key_name] = []
            
            with open(file,"r") as f:
                txt = f.read()
                
                # 以焓值作为分割线
                pattern = "((?:\n|^)-?[0-9]+\.[0-9]+\n)"
                re_lst = re.split(pattern,txt, maxsplit=20000, flags=re.S, )[1:]
#                 print(re_lst)
                
                for idx,cnt in enumerate(re_lst):
                    # POSCAR的第一行为cnt，也即焓值
                    if idx%2 == 0:
                        self.enthalpy_dict[key_name].append(float(cnt.strip()))
                    # POSCAR除了第一行的的内容为cnt
                    else:
                        # calypso中的poscar是缺少元素行的，添加上
                        pattern = "((?:\s+\d\s)+)" # eg: "2 4"
                        cnt=re.sub(pattern, f"\n   {self.element[0]}   {self.element[1]} \\1", cnt)
                        self.atoms_dict[key_name].append(cnt)
    
    def save_as_poscar(self):
        """将POSCAR文本储存为.vasp文件"""
        self.file_paths = []
        for key in self.atoms_dict:
            for idx,atoms in enumerate(self.atoms_dict[key]):
                path = f"./strus/{key}_{idx+1}_{str(self.enthalpy_dict[key][idx])}.vasp"
                self.file_paths.append(path)
                with open(path,"w") as f:
                    f.write(atoms)
    @staticmethod
    def get_sym(atoms, prec=0.1):
        """获得结构的对成性"""
        return spg.get_spacegroup(cell=atoms,symprec=prec)
    
    def file2ase(self):
        """将.vasp文件转换成ase.atoms类"""
        self.atoms_list = []
        for file in self.file_paths:
            atoms = read(file, format="vasp")
            caly_index = file.split("_")[2]  # 迭代次数，比如caly_opt_18，caly_index为18
            index = file.split("_")[3]  # 当前迭代次数中的第 index个文件
            enthalpy = float(file.split("_")[-1][:-5])
            atoms.info={}
            atoms.info["caly_index"]=caly_index
            atoms.info["index"]=index
            atoms.info["formula"] = atoms.symbols.get_chemical_formula()
            atoms.info["enthalpy"]=enthalpy
            atoms.info["spacegroup"]=File2class.get_sym(atoms, prec=self.prec)
            self.atoms_list.append(atoms)
        
        return self.atoms_list
        
class GetEhull(object):
    """获得ehull"""

    def __init__(self, atoms_list, symbols_list=["N","Si"], ground_state={"N":1,"Si":1}):
        self.symbols_list = symbols_list
        # * boundary结构的ehull值，以及对应的atoms结构，后面会更新
        # * eg: boundary_ehth = {"N":1000; "P":1000}
        #self.boundary_enth = {i: 1000 for i in self.symbols_list}
        self.boundary_enth = ground_state
        # * eg: boundary_enth_strus = {"N": <type 'atoms1'>, "P": <type 'atoms2'>}
        #self.boundary_enth_strus = {i: None for i in self.symbols_list}
        # ! 去除了P1和重复结构的列表 - 最终要的结果
        self.list_prop_strus = atoms_list

    def remove_repeated_strus(self):
        """去除重复的结构
        Parameters
        ----------
        strus_list - list of ase.atoms
            存放结构的列表

        """
        ret = []
        for stru in self.list_prop_strus:
            if stru not in ret:
                ret.append(stru)
        self.list_prop_strus = ret


    def remove_P1_strus(self):
        pop_index=[]
        for idx,stru in enumerate(self.list_prop_strus):
            if spg.get_symmetry_dataset(stru, symprec=0.1)['number'] == 1:
                pop_index.append(idx)
        self.list_prop_strus = np.delete(np.array(self.list_prop_strus),pop_index).tolist()


    def updata_boundary(self):
        """更新边界（单一组分）的结构和焓值信息，该信息用来计算ehull"""
        for atoms in self.list_prop_strus:
            symbols = atoms.get_chemical_symbols()  # ["P","N","N",...]
            if len(set(symbols)) == 1:
                symbol = symbols[0]
                if atoms.info["enthalpy"] < self.boundary_enth[symbol]:
                    self.boundary_enth[symbol] = atoms.info["enthalpy"]
                    self.boundary_enth_strus[symbol] = atoms

    def update_ehulls(self):
        """更新ehull"""
        for atoms in self.list_prop_strus:
            n_atoms = len(atoms.symbols)    # 原子数 [7,7,7,6,....]
            formula = atoms.info["formula"]  # 分子数，比如"N3P3"
            enthalpy_total = atoms.info["enthalpy"] * n_atoms
            # 获得atoms每种元素的数量 
            n_symbols = Counter(atoms.get_chemical_symbols()) # {"N":5, "P":4}

            # AxBy=xA+yB，ehull=(enthalpy(AxBy)-(x*enthalpy(A)+y*enthalpy(B)))/(x+y)
            # 其中 minuend 就是 x*enthalpy(A)+y*enthalpy(B)
            minuend = n_symbols[self.symbols_list[0]] * self.boundary_enth[self.symbols_list[0]
                                                                           ] + n_symbols[self.symbols_list[1]] * self.boundary_enth[self.symbols_list[1]]
            ehull = (enthalpy_total - minuend)/n_atoms
            atoms.info["ehull"] = ehull
    
    @staticmethod
    def print_info(info, index):
        print(f"index\t 分子式\t  空间群\t enthalpy  \tehull \t 迭代次数\t索引\t in_ehull")
        for idx,i in enumerate(info):
            print(f"{idx:<8d}{i[0]:<10s}{i[1]:<15s}{i[2]:<12f}{i[3]:<15f}{i[4]:<10s}\t{i[5]:<10s}{idx in index}")
                
    def plot_convex_hull(self, show=True):
        refs = []  # eg: [("N3",ehull1),("N3P3",ehull2)]
        info = []
        for atoms in self.list_prop_strus:
            formula = atoms.symbols.get_chemical_formula()
            ehull = atoms.info["ehull"]
            refs.append((formula, ehull))
            
            spacegroup = atoms.info["spacegroup"]
            caly_index = atoms.info["caly_index"]
            index = atoms.info["index"]
            formula = atoms.info["formula"]
            enthalpy = atoms.info["enthalpy"]
            
            info.append((formula, spacegroup, enthalpy, ehull, caly_index, index))
        
        # 在ehull上的结构索引
        self.pd = PhaseDiagram(refs,verbose=False)
        ehull_index = np.where(self.pd.hull == True)[0]   
        GetEhull.print_info(info, ehull_index)
        self.pd.plot()
        
    def save_ehull_strus_to_vasp(self):
        """保存在ehull线上的结构"""
        index = np.where(
            self.pd.hull == True)[0]  # 判断结构是否在ehull线上 eg: [[True, False, ....]]
        for i in index:
            write(f"{i}.vasp", self.list_prop_strus[i], format="vasp")

            
# 创建类
f = File2class("./",["N","Si"],0.1)
# 获得文件路径
f.calc_path()
# 获得指定文件路径
f.calc_dest_path()
# 获得POSCAR文本
f.calc_atoms_dict()
# 保存vasp文件
f.save_as_poscar()
# 将vasp文件转换成 ase.atoms 类
atoms_list = f.file2ase()

e = GetEhull(atoms_list, ["N","Si"], {"N":-6.60416,"Si":-3.82484}) # 给定基态能量
e.remove_repeated_strus()
e.remove_P1_strus()
e.updata_boundary()
e.update_ehulls()
e.plot_convex_hull()
e.save_ehull_strus_to_vasp()