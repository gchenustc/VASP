from ase.io import read
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt
import numpy as np
import sys

"""
desc: 1. calculating the average bond length along the z axis and plot the line graph. each y value represents the average appointed typed bond length less the corresponding x value(that is, structure's z length) 2. the params is the path for the structures. you can pass more than 1 structure and not more than 7.
3. the bond type you can rewrite in the following content. 4. in the current version, the script can't give the graph data.
author: gchen
time: 20240118
"""

#### setting ####
POSCARS_PATH = sys.argv[1:]
STA_ATOMIC_TYPE = "N"  # 要统计的原子类型
BOND_TYPE = ["N", STA_ATOMIC_TYPE]  # 计算的何种原子间的键长
#### setting ####

def getBondLenList(stru):
    """
    para: ase struture type
    return: type is list in the list, inner list store the bond length and inner list index is the atomic index.
    """
    n_sta = stru.get_chemical_symbols().count(STA_ATOMIC_TYPE)  # 统计的原子的数量
    # 统计的原子所连接的其他原子列表, [[1.43,1.5,1.42],[1.35, ...], ...]
    atomic_bondLen_list = [[] for i in range(n_sta)]

    ana = Analysis(stru)
    # [[(index1,index2),(..., ...), ...]]
    bondIndexList = ana.get_bonds(BOND_TYPE[0], BOND_TYPE[1], unique=True)

    bondLenthList = ana.get_values(bondIndexList)[0]  # [len1, len2, len3, ...]

    for idx, index_pair in enumerate(bondIndexList[0]):
        for index in index_pair:
            atomic_bondLen_list[index].append(bondLenthList[idx])

    return atomic_bondLen_list


def filterSelectedAtom(stru):
    symbols = stru.get_chemical_symbols()
    indices = [i for i, x in enumerate(symbols) if x == STA_ATOMIC_TYPE]
    atoms_seleted_type = stru[indices]
    return atoms_seleted_type


def getAverageBondLen(atomic_bondLen_list):
    return [sum(i)/len(i) for i in atomic_bondLen_list]


def getAverageValueIn2dlist(_2dlist):
    # 计算二维列表的平均值
    total = 0
    count = 0
    for row in _2dlist:
        for num in row:
            total += num
            count += 1

    return total / count


def getAverageBondLenAccordingZ(atomic_bondLen_list):
    atomic_average_bondLen_list = []
    for idx, _ in enumerate(atomic_average_bondLen):
        tmp = getAverageValueIn2dlist(atomic_bondLen_list[:idx+1])
        atomic_average_bondLen_list.append(tmp)
    return atomic_average_bondLen_list


if __name__ == "__main__":
    alphas = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    lw = [4, 3.5, 3, 2.5, 2, 1.5, 1]

    for idx, poscar_path in enumerate(POSCARS_PATH):
        stru_all = read(poscar_path, format="vasp")  # read structure

        # rank the atom sequence according z value.
        indices = list(map(lambda x: x.index, sorted(
            stru_all, key=lambda x: x.z)))  # 获得排序后的原子索引
        stru_all = stru_all[indices]

        # filter selected atom
        stru_selected = filterSelectedAtom(stru_all)

        atomic_bondLen_list = getBondLenList(stru_selected)
        atomic_average_bondLen = getAverageBondLen(atomic_bondLen_list)

        atomic_average_bondLen_according_z = getAverageBondLenAccordingZ(
            atomic_bondLen_list)
        # print(atomic_average_bondLen_according_z)

        x = list(map(lambda x: x.z, stru_selected))
        y = atomic_average_bondLen_according_z

        plt.plot(x, y, label=poscar_path, lw=lw[idx], alpha=alphas[idx])
        plt.xlabel("z axis distance")
        plt.ylabel("average N-N bond lenth(Int)")
        plt.legend(frameon=False)
        plt.savefig("AvebondLengthAsZaxis.png")
