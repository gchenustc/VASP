import numpy as np
import os
import argparse
from ase.io.pov import get_bondpairs
from ase.io import read
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms

"""
******** 使用方法 *********

0. 配置好ase执行环境
1. 在脚本里调整vasp参数
2. 运行脚本 
python current_file.py -l lobster  # -l 后指定lobster的执行命令
python current_file.py -l lobster -p -a 2 4 5 7  # 出现 -p 则表示计算指定原子的cohp(否则计算所有原子的cohp), -a 后为指定原子的序号，比如示例将计算序号为2和4以及5和7原子的cohp
******** END 使用方法 *********
"""

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
      directory="./workdir", # 计算路径
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
      kpts = [13,13,1],
      gamma = "True",
      nbands = 242,

      #algo = "FAST",
      ialgo=38,
      nelm=300,
      nelmin=8,
      ediff=1E-6,
      ismear=0,
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
# calc = Vasp(restart=True, directory="vasp")  # ../vasp为vasp存储路径
atoms = calc.get_atoms()

# 获得成键原子的索引
# [(1, 31, array([0, 0, 0])), ..., ]
bondpairs_raw = get_bondpairs(atoms) 

bondpairs = []
for bp in bondpairs_raw:
    tmp = list(bp[:2])  # [1,31], ...
    tmp.sort() # 排序，因为还存在[31,1]，排除等价情况
    if tmp not in bondpairs:
        bondpairs.append(tmp)
    bondpairs.sort(key=lambda x:x[0]) # 按照list里的每个子list的第一个数排序 --> [[0,30],...,[1,31],...,]

f = open('labels', 'w')
for bp in bondpairs:
    f.write(f"{atoms.get_chemical_symbols()[bp[0]]}[{bp[0]}]-{atoms.get_chemical_symbols()[bp[1]]}[{bp[1]}]\n")
f.close()

# 写入lobster 的输入参数
with open('lobsterin', 'w') as f:
    f.write('COHPstartEnergy  -30\n')
    f.write('COHPendEnergy  10\n')
    f.write('usebasisset pbeVaspFit2015\n')
    f.write('useRecommendedBasisFunctions\n')
    f.write('gaussianSmearingWidth 0.05\n')
    if not args.appoint:
        for bp in bondpairs:
            f.write(f'cohpbetween atom {bp[0]+1} and atom {bp[1]+1}\n')
    else:
        # 将输入[2 4 6 3] 转换成 [[2,4],[6,3]]
        bondpairs = [(args.atoms[i],args.atoms[i+1]) for i in range(len(args.atoms)) if not i%2]
        for bp in bondpairs:
            f.write(f'cohpbetween atom {bp[0]} and atom {bp[1]}\n')

# 执行lobster      
os.system(args.lobster)
