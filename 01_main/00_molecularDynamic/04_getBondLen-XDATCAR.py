from ase.io.vasp import read_vasp_xdatcar
from ase.io import read
from ase.geometry.analysis import Analysis
from collections import Counter
import multiprocessing
import matplotlib.pyplot as plt
from ase.data import atomic_numbers, covalent_radii
from multiprocessing import Manager
import time
"""
DESC: 计算指定键长随着步数的变化。需要适当调整共价半径(用**标注了修改位置)
AUTHOR: gchen
"""

xdatcar = read_vasp_xdatcar('XDATCAR', index=slice(-1))
poscar = read("POSCAR",format="vasp")
xdatcar.insert(0,poscar)
xdatcar = xdatcar[7000:10000] #取前10000个构型
bondType = ["N","N"] # 计算的键长类型

chemSymbolsList = xdatcar[0].get_chemical_symbols() # ["N"...,"F"}]
numCounts = Counter(chemSymbolsList) # Counter({'N': 54, 'F': 54})
elements = list(set(chemSymbolsList)) # ["N","F"]
numElements = len(numCounts) # 2

# 计算N-N共价半径
kind1_number = atomic_numbers[bondType[0]]
kind2_number = atomic_numbers[bondType[1]]
kind1_radii = covalent_radii[kind1_number] * 1.06 # **
kind2_radii = covalent_radii[kind2_number] * 1.06 # **
conv_radiis = kind1_radii + kind2_radii

def process_data(stru):
    """计算平均N-N键长"""
    ana = Analysis(stru)
    bondIndexList = ana.get_bonds(bondType[0], bondType[1], unique=True)
    bondLenthList = ana.get_values(bondIndexList)[0]
    for i in range(len(bondLenthList)-1, -1, -1,):
        if bondLenthList[i] >= conv_radiis:
            del bondLenthList[i]
    return  sum(bondLenthList) / len(bondLenthList)

def write_file(name, lst):
    # 写入文件
    with open(name, 'w') as file:
        # 逐行写入内容
        for idx,i in enumerate(lst):
            file.write(f"{str(idx+1)}\t")
            file.write(f"{str(i)}\n")
        
        file.write(f"average\t{sum(lst) / len(lst)}")

def plot(name, lst, xlabel,ylabel):
    plt.plot(range(1,len(lst)+1), lst)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(name)
    
if __name__ == "__main__":
    
    start_time = time.time()  # 记录开始时间
    
    with multiprocessing.Pool(1) as pool:
        averBlList = [i for i in pool.map(process_data, [stru for stru in xdatcar])] 
#    averBlList = [process_data(stru) for stru in xdatcar] 

    end_time = time.time()  # 记录结束时间
    print("计算完成，总共用时：", end_time - start_time, "秒")
    
    write_file('averageBondLen.txt', averBlList)
    plot("averageBondLen.png", averBlList, 'step', 'average bond lenth')
    
