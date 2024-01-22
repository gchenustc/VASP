import ase
from ase.io import read
from ase.io.pov import get_bondpairs
from ase.calculators.vasp import Vasp
from ase.build.supercells import make_supercell
import numpy as np
import time
from pymatgen.core.periodic_table import Element
# print(Element("H").atomic_radius_calculated)  # 计算原子半径
"""
描述: vasp提供CONTCAR即可

需要环境: python3
包: numpy, ase

运行方法:
```
python AveBoundLen.py
```
"""
start_time = time.perf_counter()

#calc = Vasp(restart=True, directory=".")
#atoms = calc.get_atoms()

atoms = read("CONTCAR", format="vasp")
n_atoms = len(atoms.get_atomic_numbers())  # 原子数
# >>> 扩两倍胞，因为只计算在当前周期内的键长，边界上两原子成键距离会被算长
#P = np.zeros((3, 3))
#np.fill_diagonal(P, [1, 1, 1])
#atoms = make_supercell(atoms, P)
# <<<
#print(dir(atoms))

# 获得成键电子
bondpairs = get_bondpairs(atoms)
"""
bondpairs:
[(0, 30, array([0, 0, 0])),
(1, 31, array([0, 0, 0])),....]
"""

dis_list = {}
for a0, a1, index in bondpairs:
    # 下面的判断是过滤超过一个周期的原子(为什么减1？因为Python索引从0开始)
    # if a0>n_atoms-1 or a1>n_atoms-1:
    #    continue

    dis = atoms.get_distance(a0, a1)
    
    symbol1=atoms.get_chemical_symbols()[a0]
    symbol2=atoms.get_chemical_symbols()[a1]
    _bias=0.4
    bond_len_limit = Element(symbol1).atomic_radius_calculated + Element(symbol1).atomic_radius_calculated + _bias
    #print(bond_len_limit)
    if dis > bond_len_limit:
        continue
    # print(a0,a1,dis)
    key = '-'.join(sorted([atoms[a0].symbol, atoms[a1].symbol]))
    dis_list[key] = dis_list.setdefault(key, []) + [dis]

ave = 0
n_total = sum(map(lambda x: len(x), dis_list.values()))
# print(n_total)
print("average per kind:\n")
for key, value in dis_list.items():
    mean_ = np.array(value).mean()
    #print(n_total, value)
    ave += mean_ * (len(value)/n_total)
    print(f"{key}: {mean_}")

print()

print("average for all:", ave)
time_cost = time.perf_counter() - start_time
print("---- Time used: %.4f s ----" % time_cost)
