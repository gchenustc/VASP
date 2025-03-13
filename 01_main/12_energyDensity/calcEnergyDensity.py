from pymatgen.io.vasp.outputs import Outcar,Oszicar
from pymatgen.io.vasp import Poscar
import argparse
from pymatgen.core.periodic_table import Element  
from collections import Counter
from ase.io import write,read
# 需要POSCAR和OSZICAR
# 需要的python包： pymatgen
# 用法：python 此文件


"""
比如CO2分解为CO和O2 2CO2 = 2CO+O2
则 ground_state_input = [{"C":1, "O":1},{"O":2}] 
ratio = [2,1]
energy = [1*CO的能量, 1*O2的能量] # 要和ground_state_input对应
"""
#Fibrous_P_0GPa=-5.54817
#N_alpha_0GPa=
#N_epsilon_0GPa=
#Si_0GPa=-5.424965

ground_state_input = [{"P" : 3, "N": 5},{"N": 2}] # 基态各个元素的数量  ZrN6 --> Zr + 3N2
ratio = [1,2] # 基态比例
energy_input = [-58.9965775, -16.720379]   # 基态元素的能量: Zr, N2
# 需要编辑的部分 

energy = energy_input[:]
for i in range(len(energy_input)):   
    energy[i] = energy_input[i] * ratio[i]
    
# 读取含能物质能量
oszicar=Oszicar(filename="OSZICAR") 
final_energy_product = oszicar.final_energy

# 读取结构
poscar = Poscar.from_file(filename="POSCAR")
stru = poscar.structure
stru_ase = read("POSCAR", format="vasp")

# 获得符号
symbols = stru_ase.get_chemical_symbols() # ["Zr","N","N","N","N", "N","N"]
# 统计各个原子出现次数
symbols_sta = Counter(symbols)  #Counter({'Zr': 1, 'N': 6})

from collections import Counter
# 基态统计, 经过下面计算得到该值为，Counter({'Zr': 1, 'N': 6})
ground_state_sta = Counter() 
for i,item in enumerate(ground_state_input):
    for key,value in item.items():
        item[key]=value*ratio[i] 
    ground_state_sta.update(item)

# 判定输入是否正确
keys = sorted(list(set(symbols_sta.keys()).union(set(ground_state_sta.keys()))))  # ["N", "Zr"]

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
print(f"读取的含能物质总能：{final_energy_product}")
print(f"计算的含能物质总能 (与基态物质保持相同数量) {final_energy_product}/{n} = {numerator}")
for i in range(len(ratio)):
    print(f"计算过程 (含能物质能量-基态物质能量)： {numerator}- ({energy_input[i]}) * {ratio[i]} = {numerator-energy[i]}")
    numerator -= energy[i]
    
# 分子量
mass = 0
for key in keys:
    mass +=  ground_state_sta[key] * Element(key).atomic_mass.real
print(f"含能物质的分子量: {mass}")

print(f"能量密度 = {numerator} * 96.485 / {mass} = {(numerator * 96.485) / mass}")