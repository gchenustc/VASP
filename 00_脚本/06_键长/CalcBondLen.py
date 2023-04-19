from ase.io import read
from pymatgen.core.periodic_table import Element
from ase.geometry.analysis import Analysis
from pymatgen.io.vasp import Poscar
# 计算POSCAR键长和平均键长
# 需要ase和pymatgen包
# 用法 "python 此文件" 即可
stru_ase = read("POSCAR", format="vasp")

poscar = Poscar.from_file(filename="POSCAR")
stru_pymatgen = poscar.structure

ana = Analysis(stru_ase)

_kinds = stru_pymatgen.formula.split()   #['P1', 'N7']
kinds = {x[0]: x[1:] for x in _kinds}   #{'P': '1', 'N': '7'}
n_kinds = len(kinds) # 2

bonds_index = {}
bonds_value = {}
for kind1 in kinds:
    for kind2 in kinds:
        bond_index = ana.get_bonds(kind1, kind2, unique=True)
        if not kind2+"-"+kind1 in bonds_index: # 排除相同的情况，P-N 反过来 N-P如果存在，则P-N不记录
            bonds_index[kind1+"-"+kind2] = bond_index[0]
            if bond_index[0]:  # 只有 kind1和kind2有成键的情况才可以调用get_values()
                bonds_value[kind1+"-"+kind2] = ana.get_values(bond_index)[0]
            #else:
                #bonds_value[kind1+"-"+kind2] = []
    
for key, value in bonds_value.items():
    value_str = ", ".join(str(item) for item in value)  # 用逗号分隔列表中的整数
    print(f"{key}: {value_str}")
print()
for key, value in bonds_value.items():
    average = sum(value) / len(value)  # 计算平均值
    print(f"The average of {key} is {average}")