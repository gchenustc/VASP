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
python *.py -c POSCAR1 POSCAR2 -m 0 0 1 # 将POSCAR1和POSCAR2沿着c轴移动1A，可以不传入-c，默认文件名为POSCAR
python *.py -c POSCAR1 POSCAR2 -m 0 0 1 -v # -v 可视化处理后的图像
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

if prm.visualize==True:
    view(i_poscars_r)
    
# moving operate
for position in positions:
    position += moving
    
# Atom class
atoms = [ase.Atoms(numbers=n, cell=c, positions=p, pbc=True) for n,c,p in zip(numbers, lattice, positions)]

# write
[ase.io.write(f"POSCAR_out{i+1}", atom, format='vasp') for i,atom in enumerate(atoms)]

endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print(time.strftime("%H:%M:%S on %A %d %B %Y"))
print("Program was running for %.2f seconds." % runtime)
