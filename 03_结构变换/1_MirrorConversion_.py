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
description: Mirror symmetry for correscoping axis(eg. x, y or z)
usage:
python *.py -c POSCAR1 POSCAR2 -x -y # mirror symmetry of x and y axis for cell named POSCAR1 and POSCAR2
"""

# argparser
parser = argparse.ArgumentParser(description="Mirror symmetry for correscoping axis(eg. x, y or z)")
parser.add_argument('-c','--POSCAR', default=["POSCAR"], nargs="+", help="Specify the POSCAR file")
parser.add_argument('-v','--visualize', action="store_true", default=False, help='Use ASE-GUI to visualize structure (Default=False)')
parser.add_argument('-x', action="store_true", default=False, help="mirror symmerty for x axis")
parser.add_argument('-y', action="store_true", default=False, help="mirror symmerty for y axis")
parser.add_argument('-z', action="store_true", default=False, help="mirror symmerty for z axis")
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


# mirror operate
def mirror(scaled_position, axis):
    scaled_position[:,axis] = 1 - scaled_position[:,axis]

if prm.x:
    [mirror(scaled_position, 0) for scaled_position in scaled_positions]
if prm.y:
    [mirror(scaled_position, 1) for scaled_position in scaled_positions]
if prm.z:
    [mirror(scaled_position, 2) for scaled_position in scaled_positions]

# final poscar class
atoms = [ase.Atoms(numbers=n, cell=c, scaled_positions=p, pbc=True) for n,c,p in zip(numbers, lattice, scaled_positions)]

# write
[ase.io.write("POSCAR_out", atom, format='vasp') for atom in atoms]

endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print(time.strftime("%H:%M:%S on %A %d %B %Y"))
print("Program was running for %.2f seconds." % runtime)