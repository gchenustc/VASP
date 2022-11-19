from ase.io import read
from ase.calculators.vasp import Vasp
import numpy as np
import matplotlib.pyplot as plt
# 要提前准备坐标文件
cgN = read(filename="cubic/cgN.vasp", format="vasp")

# ENCUT范围
encuts = [i for i in range(300,700,50)]

# 设置计算器
from ase.dft import kpoints as kp
kpts = kp.monkhorst_pack([10,10,10]) # 设置k点

calcs = [Vasp(directory=f"encut_{en}", # 计算路径
      txt=f"encut_{en}.log", # vasp输出文件名

      istart=0,
      icharg=2,
      xc="pbe",
      prec="Accurate",
      ediff=1E-3,
      ibrion=-1,
      ismear=0,
      sigma=0.05,
      lwave=".F.",
      lcharg=".F.",
      encut=en,  # 变量
      kpts = kpts,
      algo="very_fast") for en in encuts]

cgNs = [cgN.copy() for calc in calcs]
for index,cgN in enumerate(cgNs):
    cgN.calc = calcs[index]
    
energies = [atoms.get_potential_energy() for atoms in cgNs]

plt.plot(encuts, energies, "bo-")
plt.xlabel("ENCUT (eV)")
plt.ylabel("Total energy (eV)")
plt.savefig("cgN-encut-v.png")
