from ase.io import read
from ase.calculators.vasp import Vasp
import numpy as np
import matplotlib.pyplot as plt
# 要提前准备坐标文件
atoms = read(filename="POSCAR", format="vasp")

# ENCUT范围
encuts = list(range(200,700,20))

# 设置计算器
calcs=\
[
    Vasp(
      directory=f"encuts_{i}eV", # 计算路径
      txt="./out.log", # vasp输出文件名
      xc="pbe",
      npar=8,
      #lreal="auto",
      #isym=0,
      istart=0,
      icharg=2,
      prec="Accurate",
      ibrion=-1,
      isif=2,
      #nsw=300,
      potim=0.5,
      encut=i,
      nelm=150,
      nelmin=6,
      ediff=1E-6,
      #ediffg=-1E-3,
      ismear=0,
      sigma=0.05,
      lwave=False,
      lcharg=False,
      kspacing = 0.25,
      kgamma=False,
      algo = "FAST",
      #pstress = 2000,
      #ialgo=38,
      setups='recommended',
      ) for i in encuts
 ]

atoms_total = []
for calc in calcs:
    atoms = atoms.copy()
    atoms.calc = calc
    atoms.get_potential_energy()
    atoms_total.append(atoms)
 
energies = [atoms.get_potential_energy() for atoms in atoms_total]

plt.plot(encuts, energies, "bo-")
plt.xlabel("ENCUT (eV)")
plt.ylabel("Total energy (eV)")
plt.savefig("encut-e.png")
