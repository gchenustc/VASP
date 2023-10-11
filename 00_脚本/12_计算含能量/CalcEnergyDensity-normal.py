from pymatgen.io.vasp.outputs import Outcar,Oszicar
from pymatgen.io.vasp import Poscar
import argparse
from pymatgen.core.periodic_table import Element  
from collections import Counter
from ase.io import write,read
# 需要POSCAR和OSZICAR
# 需要的python包： pymatgen
# 用法：python 此文件


# 需要编辑的部分
"""
比如 2CO+O2=2CO2
则 ground_state_comb = [{"C":1, "O":1},{"O":2}] 
ratio = [2,1]
energy = [1*CO的能量, 1*O2的能量] # 要和ground_state_comb对应
"""
#Fibrous_P_0GPa=-5.54817
#N_alpha_0GPa=-8.37770
#N_epsilon_0GPa=-8.37539
#Si_0GPa=-5.424965

ground_state_comb = [{"N": 1},{"Si":1}] # 基态各个元素的数量
ratio = [6,1] # 基态比例
energy=[-8.37539,-5.424965]   # 基态元素的能量
# 需要编辑的部分 


for i in range(len(energy)):
    energy[i] = energy[i] * ratio[i]
# 读取能量
oszicar=Oszicar(filename="OSZICAR")
final_energy_product = oszicar.final_energy

# 读取结构
poscar = Poscar.from_file(filename="POSCAR")
stru = poscar.structure
stru_ase = read("POSCAR", format="vasp")

# 获得符号
symbols = stru_ase.get_chemical_symbols() # ["F","N","N",...]
# 统计各个符号出现次数
symbols_sta = Counter(symbols)  #Counter({'F': 8, 'N': 8})

from collections import Counter
ground_state_sta = Counter() # 基态统计  Counter({'N': 3, 'F': 3})
for i,item in enumerate(ground_state_comb):
    for key,value in item.items():
        item[key]=value*ratio[i] 
    ground_state_sta.update(item)
    
# 判定输入是否正确
keys = sorted(list(set(symbols_sta.keys()).union(set(ground_state_sta.keys())))) #["F,"N"]
ratio_ = symbols_sta[keys[0]] / ground_state_sta[keys[0]]
for key in keys:
    if symbols_sta[key] / ground_state_sta[key] != ratio_:
        print("基态比例设置错误")
        raise KeyError
        
        
# 确定末态物质数量
n = symbols_sta[keys[0]] / ground_state_sta[keys[0]]

# 计算含能量
# 分子
numerator = final_energy_product/n
for i in range(len(ratio)):
    numerator -= energy[i]
    
# 分子量
mass = 0
for key in keys:
    mass +=  ground_state_sta[key] * Element(key).atomic_mass.real

print("energy density:", numerator*96.485/mass)