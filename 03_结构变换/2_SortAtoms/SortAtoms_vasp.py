from ordered_set import OrderedSet
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
description: Sort atoms along correscoping axis(eg. x, y or z)
usage:
python *.py -c POSCAR1 POSCAR2 -x
"""

# argparser
parser = argparse.ArgumentParser(description="Sort atoms along the specific axis")
parser.add_argument('-c','--POSCAR', default=["POSCAR"], nargs="+", help="Specify the POSCAR file")
parser.add_argument('-v','--visualize', action="store_true", default=False, help='Use ASE-GUI to visualize structure (Default=False)')
group = parser.add_mutually_exclusive_group()
group.add_argument('-x', action="store_true", default=False, help="x axis sorting")
group.add_argument('-y', action="store_true", default=False, help="y axis sorting")
group.add_argument('-z', action="store_true", default=False, help="z axis sorting")
prm = parser.parse_args()

# Starting
starttime = time.time()
print("Starting calculation at", end='')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))

# strip
i_poscars = [pos.strip() for pos in prm.POSCAR]
i_poscars_r = [ase.io.read(pos) for pos in i_poscars]

# get cell data
n = [pos.get_global_number_of_atoms() for pos in i_poscars_r] # 原子数
numbers = [pos.numbers for pos in i_poscars_r]
positions = [pos.positions for pos in i_poscars_r]
scaled_positions = [pos.get_scaled_positions() for pos in i_poscars_r]
lattice = [pos.cell for pos in i_poscars_r]
cell_lenth_angle = [pos.get_cell_lengths_and_angles() for pos in i_poscars_r]
symbols = [pos.get_chemical_symbols() for  pos in i_poscars_r]

if prm.visualize==True:
    view(i_poscars_r)


# 排序操作
def sortAtoms(numbers, positions, axis):
    
    # 相同元素簇放一起
    indicx_temp = np.array(numbers).argsort()[::-1] # 按分子量从大到小
    numbers = numbers[indicx_temp]
    positions = positions[indicx_temp]
#     print(positions)
    
    # 对不同的簇进行排序
    groups = OrderedSet(numbers)
    n_groups = len(groups)
    n_each_groups = []
    for group in groups:
         n_each_groups.append(numbers.tolist().count(group)) 
#     print(n_each_groups) # [3,1]
    
    index = n_each_groups
    index.insert(0,0)
    index = np.array(index).cumsum()
#     print(index) # [0,3,4]
    for i in range(1,len(index)):
        i_section = positions[index[i-1]:index[i]]
        f_section = i_section[np.argsort(i_section[:,axis])]
        positions[index[i-1]:index[i]]= f_section
    return numbers, positions

if prm.x:
    axis = 0
elif prm.y:
    axis = 1
else:
    axis = 2
numbers,positions = zip(*[sortAtoms(n,p,axis) for n,p in zip(numbers, positions)])

# final poscar class
atoms = [ase.Atoms(numbers=n, cell=c, positions=p, pbc=True) for n,c,p in zip(numbers, lattice, positions)]

# write
[ase.io.write("POSCAR_out", atom, format='vasp') for atom in atoms]

endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print(time.strftime("%H:%M:%S on %A %d %B %Y"))
print("Program was running for %.2f seconds." % runtime)