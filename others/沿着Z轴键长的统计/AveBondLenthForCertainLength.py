from ase.io import read
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt
import numpy as np
import sys
"""
desc: 1. calculating the average bond length along the z axis and plot the scatter with line graph. each y value represents the average appointed typed bond length for the centain range which less than x value (that is, structure's z length). the range depends the division number. For example. if the z length equals 10 A, and the NUM_DIVISION equals 5, each y value is the average bond-length for 2A's distances less than correspongding x value. 2. the params is the path for the structures. you can pass more than 1 structure and not more than 7.
3. the bond type you can rewrite in the following content. 4. in the current version, the script can't give the graph data.
author: gchen
time: 20240118
"""

#### setting ####
POSCARS_PATH = sys.argv[1:]
STA_ATOMIC_TYPE = "N"  # 要统计的原子类型
BOND_TYPE = ["N", STA_ATOMIC_TYPE]  # 和谁的键长
NUM_DIVISION = 10 # x轴点的数量
#### setting ####

def filterSelectedAtom(stru):
    symbols = stru.get_chemical_symbols()
    indices = [i for i, x in enumerate(symbols) if x == STA_ATOMIC_TYPE]
    atoms_seleted_type = stru[indices]
    return atoms_seleted_type

def getAverageBondLen(atomic_bondLen_list):
    return [sum(i)/len(i) for i in atomic_bondLen_list]

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


def getAverageValueIn2dlist(_2dlist):
    # 计算二维列表的平均值
    total = 0
    count = 0
    for row in _2dlist:
        for num in row:
            total += num
            count += 1

    return total / count

def getValueIndexInCoordsList(value, coords):
    """
    获得value在列表coords中的索引
    coords中的元素是二维坐标，如果value在范围内，则返回该范围所在的index
    """
    for idx,coord in enumerate(coords):
        if coord[0] <= value<= coord[1]:
            return idx
    return None

def get_indexes(lst, element):
    indexes = [i for i, x in enumerate(lst) if x == element]
    return indexes

def getAverageBondLenForCertainLen(stru_selected, atomic_bondLen_list):
    distance_furthest = stru_selected[-1].z
    endpoint = [(distance_furthest/NUM_DIVISION) * i for i in range(NUM_DIVISION+1)]
    """
    endpoint :[0.0, 2.4604009750697586, 4.920801950139517, 7.381202925209276, 9.841603900279035, 12.302004875348793]
    """
    endpoint_pair = []
    """
    endpoint_pair :[(0.0, 2.4604009750697586), (2.4604009750697586, 4.920801950139517), (4.920801950139517, 7.381202925209276), (7.381202925209276, 9.841603900279035), (9.841603900279035, 12.302004875348793)]
    """
    for i in range(len(endpoint)-1):
        tmp = endpoint[i],endpoint[i+1]
        endpoint_pair.append(tmp)
    
    #centerpoint = [(i[0]+i[1])/2 for i in endpoint_pair]
    
    distance_z_list = list(map(lambda x: x.z, stru_selected))
    atomic_division_index = []
    
    for z in distance_z_list:
        index = getValueIndexInCoordsList(z, endpoint_pair)
        atomic_division_index.append(index)

    ave_value = []
    for i,_ in enumerate(endpoint[1:]):
        indexes = get_indexes(atomic_division_index, i)
        tmp = np.array(atomic_bondLen_list,dtype=object)[indexes].tolist() # 由索引获得value
        ave_value.append(getAverageValueIn2dlist(tmp))
    
    return endpoint[1:], ave_value
    
def write_file(name, x, y):
    # 写入文件
    with open(name, 'w') as file:
        # 逐行写入内容
        for i,j in zip(x,y):
            file.write(f"{str(i)}\t")
            file.write(f"{str(j)}\n")
            
if __name__ == "__main__":
#    alphas = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
#    lw = [4, 3.5, 3, 2.5, 2, 1.5, 1]
#    markers = ['s', 'o', 'H', '*', '^', 'v']   # 不同的样式
#    sizes = [10, 8.5, 7, 5.5, 4, 3]  # 不同的大小

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

        x,y = getAverageBondLenForCertainLen(stru_selected, atomic_bondLen_list)
        
        write_file(f"AveBondLenthForCertainLength-{poscar_path}.txt", x, y)

        #plt.plot(x, y, label=poscar_path, marker=markers[idx], markersize=sizes[idx], linestyle='-', lw=lw[idx], alpha=alphas[idx])
        plt.plot(x, y, label=poscar_path, marker="o", markersize=7, linestyle='-', lw=2)
        plt.xlabel("z axis distance")
        plt.ylabel("average N-N bond lenth(Int)")
        plt.legend(frameon=False)
        plt.savefig("AveBondLenthForCertainLength.png")
        
    
