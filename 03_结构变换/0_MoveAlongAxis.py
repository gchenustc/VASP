import argparse
import spglib
import ase
from ase.io import read, write
from ase import Atoms
from ase.visualize import view
import numpy as np
import os
import sys
import datetime
import time
import numpy as np
"""
description: Move atoms along the specific axis
usage:
python *.py -c POSCAR1 POSCAR2 -m 0 0 1 # 将POSCAR1和POSCAR2沿着z轴（Cartesian坐标）移动1A，可以不传入-c，默认文件名为POSCAR
python *.py -c POSCAR1 POSCAR2 -d direct -m 0 0 1 -v # -d指移动mode，沿着轴移动，而不是沿着固定的笛卡尔坐标移动，如果是胞是四方结构，开不开此项都一样。
python *.py -c POSCAR1 POSCAR2 -d direct -p -m 0 0 0.1 -v # -p打开后，按照胞长的比例移动，-p参数需要-d direct
"""

def is_number(s):
    """如果字符串s为浮点数或者整数，返回True"""
    if str(s).isnumeric(): # isdigit(), isdecimal()
        return True
    try:
        float(s)
        return True
    except:
        pass
    return False

# argparser
parser = argparse.ArgumentParser(description="Move atoms along the specific axis")
parser.add_argument('-c','--POSCAR', default=["POSCAR"], nargs="+", help="Specify the POSCAR file")
parser.add_argument('-v','--visualize', action="store_true", default=False, help='Use ASE-GUI to visualize structure (Default=False)')
parser.add_argument('-d','--mode', action="store", default="cartesian", choices=["c","C","D","d","direct","cartesian","Direct","Cartesian"], help='move mode')
parser.add_argument('-p','--prop', action="store_true", help='switch moved by percentage')
parser.add_argument('-m','--moving', required=True, nargs=3, help="Specify the length moved along corresponding basis vetor(x,y,z)")
prm = parser.parse_args()

# Starting
starttime = time.time()
print("Starting calculation at", end='')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))

# 检测输入
for i in prm.moving:
    if not is_number(i):
        print("arg following -m must be the number")
moving = list(map(lambda x:float(x), prm.moving))

# strip
i_poscars = [pos.strip() for pos in prm.POSCAR]
i_poscars_r = [ase.io.read(pos) for pos in i_poscars]

n = [pos.get_global_number_of_atoms() for pos in i_poscars_r] # 原子数
numbers = [pos.numbers for pos in i_poscars_r]
positions = [pos.positions for pos in i_poscars_r]
lattice = [pos.cell for pos in i_poscars_r]
symbols = [pos.get_chemical_symbols() for  pos in i_poscars_r]
# 胞长
cell_len = [np.linalg.norm(pos.get_cell(), axis=0) for pos in i_poscars_r]

if prm.visualize==True:
    view(i_poscars_r)

# moving operate
if prm.mode[0].lower() == "d":
    if not prm.prop: # prm.moving 不是比例而是长度时
        #print(1)
        moving = np.array(moving)/np.array(cell_len)
    #print(type(i_poscars_r))
    #print(i_poscars_r[0].get_scaled_positions())
    #print(cell_len)
    #print(moving)
    #i_poscars_r[0].set_scaled_positions(i_poscars_r[0].get_scaled_positions()+moving)
    [pos.set_scaled_positions(pos.get_scaled_positions()+moving) for index,pos in enumerate(i_poscars_r)]
    #print(i_poscars_r[0].get_scaled_positions())

elif prm.mode[0].lower() == "c": 
    for position in positions:
        position += moving

else:
    print("input error")
    raise InterruptedError
    
# Atom class
atoms = [ase.Atoms(numbers=n, cell=c, positions=p, pbc=True) for n,c,p in zip(numbers, lattice, positions)]

# write
[ase.io.write(f"POSCAR_out{i+1}", atom, format='vasp') for i,atom in enumerate(atoms)]

endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print(time.strftime("%H:%M:%S on %A %d %B %Y"))
print("Program was running for %.2f seconds." % runtime)
