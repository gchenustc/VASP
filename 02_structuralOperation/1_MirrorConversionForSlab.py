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
description: Mirror symmetry for z axis of slab structure. The vacuum direction must be in the z axis.
usage:
python *.py -c POSCAR # mirror symmetry of z axis.
"""

# argparser
parser = argparse.ArgumentParser(description="Mirror symmetry for correscoping axis(eg. x, y or z)")
parser.add_argument('-c','--POSCAR', default="POSCAR", help="Specify the POSCAR file")
parser.add_argument('-v','--visualize', action="store_true", default=False, help='Use ASE-GUI to visualize structure (Default=False)')
prm = parser.parse_args()

# Starting
starttime = time.time()
print("Starting calculation at", end='')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))

# strip
i_poscar = prm.POSCAR
i_poscar_r = ase.io.read(i_poscar)

# get cell data
n = i_poscar_r.get_global_number_of_atoms()
numbers = i_poscar_r.numbers
positions = i_poscar_r.positions 
scaled_positions = i_poscar_r.get_scaled_positions()
lattice = i_poscar_r.cell
cell_lenth_angle = i_poscar_r.get_cell_lengths_and_angles() 
symbols = i_poscar_r.get_chemical_symbols()

max_z_pos = np.amax(positions, axis=0)[2]
min_z_pos = np.amin(positions, axis=0)[2]
center_pos = (max_z_pos - 0)/2
print(center_pos)

if prm.visualize==True:
    view(i_poscars_r)


# mirror operate
def mirror(positions):
    for pos in positions:
        ori_pos = pos[2]
        pos[2] = ori_pos - 2 * (ori_pos - center_pos) + 0.01

mirror(positions)

# final poscar class
atoms = ase.Atoms(numbers=numbers, cell=lattice, positions=positions, pbc=True) 


# write
ase.io.write("POSCAR_out", atoms, format='vasp')

endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print(time.strftime("%H:%M:%S on %A %d %B %Y"))
print("Program was running for %.2f seconds." % runtime)