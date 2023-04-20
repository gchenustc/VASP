from pymatgen.io.vasp.outputs import Outcar,Oszicar
from pymatgen.io.vasp import Poscar
import argparse
from pymatgen.core.periodic_table import Element  
# 需要POSCAR和OSZICAR
# 需要的python包： pymatgen
# 用法：python 此文件

# 需要编辑的部分 --> 给定计算形成能所需要反应物单质的能量-要平均到每个原子
Fibrous_P_0GPa=-5.54817
N_alpha_0GPa=-8.37770
N_epsilon_0GPa=-8.37539
energy_reaction={"N": N_alpha_0GPa, "P":Fibrous_P_0GPa}
# 需要编辑的部分 

# 定义Parse
#parser = argparse.ArgumentParser(description="None")
#parser.add_argument('-p','--poscar', default="POSCAR", help="Specify the POSCAR file")
#prm = parser.parse_args()

# 根据OUTCAR读取能量
oszicar=Oszicar(filename="OSZICAR")
final_energy_product = oszicar.final_energy
#print(final_energy_product)

# 读取POSCAR结构
poscar = Poscar.from_file(filename="POSCAR")
stru = poscar.structure

# 分析POSCAR元素类型和个数
_kinds = stru.formula.split()   #['P1', 'N7']
kinds = {x[0]: x[1:] for x in _kinds}   #{'P': '1', 'N': '7'}
n_kinds = len(kinds) # 2

# 计算含能量
def energy_content(final_energy_product, kinds, energy_reaction):
    molecular_mass = 0  # 分子量
    numerator = final_energy_product # 分子
    for key in kinds:
        molecular_mass += Element(key).atomic_mass.real * int(kinds[key])
        numerator -= int(kinds[key]) * energy_reaction[key]  # 分子
    #print(numerator)
    #print(f"生成物能量：{final_energy_product}; 分子量：{molecular_mass}; 形成能：{numerator}")
    return numerator * 96.485 / molecular_mass

# 打印含能量
print(energy_content(final_energy_product, kinds, energy_reaction))
