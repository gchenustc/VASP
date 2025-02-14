import ase
from ase import io
import argparse
import numpy as np
import sys
# Author: gchen
# Date: 20250103
# Purpose: To generate hydrostatic pressure for strus

# ----------------- 1. setting -----------------
in_stru_path = '../POSCAR'  # structure path
out_stru_path = 'POSCAR'  # output structure file
direction = 0      # 0, 1, 2 for x, y, z
scale = float(sys.argv[1])       # scale for the strain (ratio, -0.05 means reduce 5% of specified direction)
pert_numb = 1      # number of perturbation
pert_atom = 0.2    # maximum perturbation for atom (Angstrom)
pert_cell = 0.01   # maximum perturbation for cell (ratio)
# ----------------- 1. setting -----------------

def gen_per_stru(idx=0):
    # Read in the structure
    atoms = io.read(in_stru_path, format='vasp')

    # Get the lattice
    lattice = atoms.get_cell()
    lattice_dest_dir = lattice[direction]

    #print(lattice_dest_dir)

    # change the lattice
    lattice_dest_dir = lattice_dest_dir * (1 + scale)
    lattice[direction] = lattice_dest_dir

    # change the lattice
    for i in range(3):
        for j in range(3):
            lattice[i][j] = lattice[i][j] * (1 + pert_cell * np.random.uniform(-1, 1))
    
    # 在设置胞之前获得原子的分数坐标:
    positions=[]
    for i in range(len(atoms)):
        positions.append(atoms[i].scaled_position)
    
    # set the lattice
    atoms.set_cell(lattice)
    #print(atoms.get_cell())
    
    # set the positions
    atoms.set_scaled_positions(positions)

    # change the atom position
    for i in range(len(atoms)):
        for j in range(3):
            # print(atoms[i].position[j])
            atoms[i].position[j] = atoms[i].position[j] + pert_atom * np.random.uniform(-1, 1)
            # print(atoms[i].position[j])
            # print('---')

    io.write(out_stru_path, atoms, format='vasp', direct=True)

def main():
    for i in range(pert_numb):
        gen_per_stru(idx=i)


if __name__ == '__main__':
    main()
