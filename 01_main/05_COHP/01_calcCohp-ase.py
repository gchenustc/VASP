"""
Desc:
    This script allows you to calculate the COHP (Crystal Orbital Hamilton Population) using VASP (Vienna Ab initio Simulation Package) and analyze the results using the Lobster tool. COHP provides insights into the bonding and interactions between atoms in a crystal structure.
Usage:
    To use this script, follow the steps below:
    Set up the ASE (Atomic Simulation Environment) execution environment to ensure compatibility.
    Adjust the VASP parameters within the script according to your specific simulation setup.
    Run the script using the following command:
    python current_file.py -l lobster: Calculates the COHP for all atoms using the specified Lobster execution command.
    python current_file.py -l lobster -p -a 2 4 5 7: Calculates the COHP for specific atoms using the specified Lobster execution command. The -p flag indicates that only specified atoms will be considered, and the -a flag is followed by a list of atom indices to calculate the COHP for.
    After running the script, it will perform the COHP calculations using VASP and Lobster based on the provided parameters.
    The results will be available for further analysis and interpretation.
"""

import numpy as np
import os
import argparse
from ase.io.pov import get_bondpairs
from ase.io import read
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms


parser = argparse.ArgumentParser(description="None")
parser.add_argument("-l","--lobster", default="lobster", help="lobster implement command")
parser.add_argument("-p","--appoint", action="store_true", help="appoint atompairs to calculate cohp")
parser.add_argument("-a","--atoms", type=int, nargs="*", help="args \"-p\" need open, after appoint atom indexs, please notice the inputs need the a multiple of 2")
args = parser.parse_args()

if args.appoint and len(args.atoms) % 2:
    print("the args:'-a' needs appear in pairs")
    exit()

# 创建结构
atoms = read(filename="POSCAR", format="vasp")

# z轴长度
z_long = np.linalg.norm(atoms.cell[2])

# 固定原子
fixed = FixAtoms(indices=[atom.index for atom in atoms if atom.c*z_long<7])   #atom的z轴长度小于<7即为所需固定的原子,atom.c获取的是分数长度，所以要乘以胞长
atoms.set_constraint(fixed)

# 调用VASP
calc = Vasp(\
      directory="./", # 计算路径
      txt="./out.log", # vasp输出文件名

      npar=8,
      lreal="auto",
      istart=0,
      icharg=2,

      ibrion=-1,
      isif=0,
      #nsw=500,
      #ediffg=ediffg,

      xc="pbe",
      isym=-1,
      prec="Accurate",
      encut=520,
      kpts = [10,10,1],
      gamma = "False",
      nbands = 250,

      #algo = "FAST",
      ialgo=38,
      nelm=300,
      nelmin=8,
      ediff=1E-6,
      ismear=-5,
      sigma=0.05,

      lwave=".T.",
      lcharg=".F.",
      setups='recommended',\
      )

# 把结构和VASP关联
atoms.calc = calc

# run vasp
atoms.get_potential_energy()

# 用ase加载已经完成计算的vasp
# calc = Vasp(restart=True, directory="./")  # ../vasp为vasp存储路径
atoms = calc.get_atoms()

# 获得成键原子的索引
# [(1, 31, array([0, 0, 0])), ..., ]
bondpairs_raw = get_bondpairs(atoms) 

if args.appoint:
    # 将输入[2,4,6,3,...] 转换成 [[1,3],[5,2],...]
    bondpairs = [(args.atoms[i]-1,args.atoms[i+1]-1) for i in range(len(args.atoms)) if not i%2]
else:
    bondpairs = []
    for bp in bondpairs_raw:
        tmp = sorted(bp[:2])
        if tmp not in bondpairs:
            bondpairs.append(tmp)
        bondpairs.sort(key=lambda x:x[0]) # 按照list里的每个子list的第一个数排序 --> [[0,30],...,[1,31],...,]

with open('labels', 'w') as f:
    for bp in bondpairs:
        f.write(f"{atoms.get_chemical_symbols()[bp[0]]}[{bp[0]+1}]-{atoms.get_chemical_symbols()[bp[1]]}[{bp[1]+1}]\n")
# 写入lobster 的输入参数
with open('lobsterin', 'w') as f:
    f.write('COHPstartEnergy  -30\n')
    f.write('COHPendEnergy  10\n')
    f.write('usebasisset pbeVaspFit2015\n')
    f.write('useRecommendedBasisFunctions\n')
    f.write('gaussianSmearingWidth 0.05\n')
    for bp in bondpairs:
        f.write(f'cohpbetween atom {bp[0]+1} and atom {bp[1]+1}\n')

# 执行lobster      
os.system(args.lobster)
