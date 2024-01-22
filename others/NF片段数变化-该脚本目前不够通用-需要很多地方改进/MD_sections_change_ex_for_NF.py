from ase.io.vasp import read_vasp_xdatcar
from ase.io import read
from ase.neighborlist import NeighborList
import numpy as np
from ase.data import atomic_numbers,covalent_radii
from collections import Counter
import csv
"""
计算和判断分子动力学过程NF(NF,N2F2,N3F3....)片段数目随步数的变化
该脚本需要改进的优化的地方很多，如果需要计算其他结构片段数变化要读懂和修改整个脚本
"""


# 导入结构文件
#xdatcar=read("NF.xyz",index=":")
#xdatcar = read_vasp_xdatcar('XDATCAR', index=slice(-1))  # slice(-1) 提取所有; slice(1:10000) 提取1-9999个结构
xdatcar=read("NF.lammpstrj",index=":",format="lammps-dump-text")

# lammpstrj文件可能无法给出元素类型，假设该结构是NF，则可能在该文件中被标记为H和He
for i in xdatcar:
    for j in i:
        if j.symbol=="H":
            j.symbol="N"
        elif j.symbol=="He":
            j.symbol="F"


# 获得原子共价半径
atomic_numbers = xdatcar[0].get_atomic_numbers()
atomic_symbols = xdatcar[0].get_chemical_symbols()
radii = covalent_radii[atomic_numbers]

# 小于原子间共价半径的25％算一个片段
radii = list(map(lambda x:x+x*0.25, radii))


def _comb(lists):
    """将ase判断的相同片段结合在一起。比如[1,2,4],[1,2,5] == > [1,2,4,5]，124是同一个片段，125是一个片段，则1245是一个片段"""
    merged_lists = []
    for sublist in lists:
        merged = False
        for merged_list in merged_lists:
            if any(x in sublist for x in merged_list):
                merged_list.extend(sublist)
                merged = True
                break
        if not merged:
            merged_lists.append(sublist)
    return merged_lists


def get_neighborlist(stru, num):
    nl = NeighborList(radii,skin=0,self_interaction=True)
    nl.update(stru)

    indices_list=[]
    for i in range(0,num):
        indices, offsets= nl.get_neighbors(i)
        indices_list.append(indices)
    #     print(indices)
    #     for j, offset in zip(indices, offsets):
    #         print(atoms.positions[j] + np.dot(offset, xdatcar[0].get_cell()))



    indices_list = [list(set(lst)) for lst in indices_list]


    for i in range(5): # 重复合并多次
        indices_list=_comb(indices_list)
        indices_list = [list(set(lst)) for lst in indices_list] # 去重
    
    return indices_list

def judge_fragment(stru,lst):
    symbols = []
    for i in lst:
        symbols.append(stru[i].symbol)
    sta = Counter(symbols)
    return sta


n_N2F2_section_list = []
n_N2F2_atoms_list = []
n_N3F3_section_list = []
n_N3F3_atoms_list = []
n_N4F4_section_list = []
n_N4F4_atoms_list = []
n_N5F5_section_list = []
n_N5F5_atoms_list = []
n_N6F6_section_list = []
n_N6F6_atoms_list = []
n_chain_section_list = []
n_chain_atoms_list = []

n_N2_section_list=[]
n_N2_atoms_list=[]
n_F2_section_list=[]
n_F2_atoms_list=[]
n_NF2_section_list=[]
n_NF2_atoms_list=[]
n_N_section_list=[]
n_N_atoms_list=[]
n_F_section_list=[]
n_F_atoms_list=[]

for stru in xdatcar:
    indices_list = get_neighborlist(stru, num=len(atomic_numbers))

    n_N2F2_section=0
    n_N2F2_atoms=0
    n_N3F3_section=0
    n_N3F3_atoms=0
    n_N4F4_section=0
    n_N4F4_atoms=0
    n_N5F5_section=0
    n_N5F5_atoms=0
    n_N6F6_section=0
    n_N6F6_atoms=0
    n_chain_section=0
    n_chain_atoms=0
    
    n_N2_section=0
    n_N2_atoms=0
    n_F2_section=0
    n_F2_atoms=0
    n_NF2_section=0
    n_NF2_atoms=0
    n_N_section=0
    n_N_atoms=0
    n_F_section=0
    n_F_atoms=0
    for i in indices_list:
        sta = judge_fragment(xdatcar[0],i)
        if sta["N"]==2:
            n_N2F2_section += 1
            n_N2F2_atoms += 4
        if sta["N"]==3:
            n_N3F3_section += 1
            n_N3F3_atoms += 6
        if sta["N"]==4:
            n_N4F4_section += 1
            n_N4F4_atoms += 8
        if sta["N"]==5:
            n_N5F5_section += 1
            n_N5F5_atoms += 10
        if sta["N"]==6:
            n_N6F6_section += 1
            n_N6F6_atoms += 12
        if sta["N"]>=7:
            n_chain_section += 1 
            n_chain_atoms += sta["N"]*2
            
        if sta["N"]==2 and len(sta)==1:
            n_N2_section += 1
            n_N2_atoms += 2
        if sta["F"]==2 and len(sta)==1:
            n_F2_section += 1
            n_F2_atoms += 2
        if sta["N"]==1 and sta["F"]==2:
            n_NF2_section += 1
            n_NF2_atoms += 3
        if sta["N"]==1 and len(sta)==1:
            n_N_section += 1
            n_N_atoms += 1
        if sta["F"]==1 and len(sta)==1:
            n_F_section += 1
            n_F_atoms += 1
#    print(n_N2F2_atoms,n_N3F3_atoms,n_N4F4_atoms,n_N5F5_atoms,n_N6F6_atoms,n_chain_atoms)
    n_N2F2_section_list.append(n_N2F2_section)
    n_N2F2_atoms_list.append(n_N2F2_atoms)
    n_N3F3_section_list.append(n_N3F3_section)
    n_N3F3_atoms_list.append(n_N3F3_atoms)
    n_N4F4_section_list.append(n_N4F4_section)
    n_N4F4_atoms_list.append(n_N4F4_atoms)
    n_N5F5_section_list.append(n_N5F5_section)
    n_N5F5_atoms_list.append(n_N5F5_atoms)
    n_N6F6_section_list.append(n_N6F6_section)
    n_N6F6_atoms_list.append(n_N6F6_atoms)
    n_chain_section_list.append(n_chain_section)
    n_chain_atoms_list.append(n_chain_atoms)
    
    n_N2_section_list.append(n_N2_section)
    n_N2_atoms_list.append(n_N2_atoms)
    n_F2_section_list.append(n_F2_section)
    n_F2_atoms_list.append(n_F2_atoms)
    n_NF2_section_list.append(n_NF2_section)
    n_NF2_atoms_list.append(n_NF2_atoms)
    n_N_section_list.append(n_N_section)
    n_N_atoms_list.append(n_N_atoms)
    n_F_section_list.append(n_F_section)
    n_F_atoms_list.append(n_F_atoms)

step=1 # 步长
time = map(lambda x:x*step, range(1,len(xdatcar)+1))
with open('data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Time(fs)','N2F2', 'N4F4', 'N6F6', "chain",'N2','F2','NF2','N','F','N2F2', 'N4F4', 'N6F6', "chain",'N2','F2','NF2','N','F'])  # 写入表头
    writer.writerows(zip(time, n_N2F2_section_list, n_N4F4_section_list, n_N6F6_section_list, n_chain_section_list,n_N2_section_list,n_F2_section_list,n_NF2_section_list,n_N_section_list,n_F_section_list,\
                         n_N2F2_atoms_list, n_N4F4_atoms_list, n_N6F6_atoms_list, n_chain_atoms_list,n_N2_atoms_list,n_F2_atoms_list,n_NF2_atoms_list,n_N_atoms_list,n_F_atoms_list))  # 写入数据行