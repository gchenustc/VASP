import ase
from ase import Atoms, Atom
from ase.io import read, write
from ase.visualize import view
import time
import numpy as np
import argparse

"""
Author: gchen
Date: 20240516
Description: Generate random gas atoms in the specified region of the slab structure.
Usage:
    optional arguments:
    -h, --help            show this help message and exit
    -c POSCAR, --POSCAR POSCAR
                            Specify the POSCAR file
    -v, --visualize       Use ASE-GUI to visualize structure (Default=False)
    -s SPACE SPACE, --space SPACE SPACE
                            Specify the space for filling atoms. eg. 18 32 means filling atoms from 18 to 32 of z-axis.
    -m, --is_molecule     Specify if the gas is the molecule
    -b BOND_LENGTH, --bond_length BOND_LENGTH
                            Specify the bond length of molecule if the --is_molecule is True
    -f [FILLED_ELEMENT [FILLED_ELEMENT ...]], --filled_element [FILLED_ELEMENT [FILLED_ELEMENT ...]]
                            Specify the filled element
    -d FILLED_DENSITY, --filled_density FILLED_DENSITY
                            Specify the density of filled atoms. eg. the pure Ar:1.394 g/cm^3
    -e MIN_DISTANCE, --min_distance MIN_DISTANCE
                            Specify the minimum distance between atoms/molecules
    example:
    python path/for/script.py -c POSCAR -s 10 20 -m -b 1.3 -f H H -d 0.1 -e 2
"""

# argparser
parser = argparse.ArgumentParser(description="Generate random gas atoms in the specified region of the structure.")
parser.add_argument('-c','--POSCAR', default="POSCAR", help="Specify the POSCAR file")
parser.add_argument('-v','--visualize', action="store_true", default=False, help='Use ASE-GUI to visualize structure (Default=False)')

parser.add_argument('-s', '--space', nargs=2, help="Specify the space for filling atoms. eg. 18 32 means filling atoms from 18 to 32 of z-axis.")
parser.add_argument('-m', '--is_molecule', action="store_true", default=False, help="Specify if the gas is the molecule")
parser.add_argument('-b', '--bond_length', type=float, default=2, help="Specify the bond length of molecule if the --is_molecule is True")
parser.add_argument('-f', '--filled_element', nargs="*", default=["Ar"], help="Specify the filled element")
parser.add_argument('-d', '--filled_density', type=float, default=1.394, help="Specify the density of filled atoms. eg. the pure Ar:1.394 g/cm^3")
parser.add_argument('-e', '--min_distance', type=float, default=2.0, help="Specify the minimum distance between atoms/molecules")

prm = parser.parse_args()

# read structure
stru_name = prm.POSCAR
stru = read(stru_name, format="vasp")

# cell data
cell = stru.cell
cell_length = np.linalg.norm(cell, axis=1)

# space for filling atoms
if not prm.space:
    start = stru.get_scaled_positions().max(axis=0)[2] * cell_length[2]
    end = cell_length[2]
else:
    start = float(prm.space[0])
    end = float(prm.space[1])


is_molecule = prm.is_molecule
filled_element = prm.filled_element  # if is_molecule equals false, filled_element only have one element
bond_length = prm.bond_length # if not is_molecule, bond_length is not useful

filled_density = prm.filled_density # g/cm^3
min_distance = prm.min_distance # minimum distance between atoms/molecules


def generate_random_point(start, end):
    """随机生成一个在指定区域内的点"""
    return np.random.uniform(low=[0,0,start], high=[1,1,end])

def is_appropriate_pos(filled_atoms, stru, min_distance, is_molecule = False):
    """检查新点(索引为-1的原子)是否位置合适"""
    
    temp = Atoms(cell=stru.cell, pbc=True)
    temp.extend(stru)
    temp.extend(filled_atoms)
    
    n = len(temp)
    
    if is_molecule:
        # if there are the second atom in molecule, don not test the neighbour atom.
        dist_list = temp.get_distances(n-1, list(range(n-2)), mic=True)
    else:
        dist_list = temp.get_distances(n-1, list(range(n-1)), mic=True)
    
    # print(dist_list)
    
    if any(np.array(dist_list) < min_distance):
        return False
    
    return True


def main():
    if is_molecule and len(filled_element) != 2:
        print("Error: The molecule should contain two atoms.")
        return False

    if (not is_molecule) and len(filled_element) != 1:
        print("Error: The filled_element should be a string.")
        return False

    filled_atom_number = ase.data.atomic_numbers[filled_element[0]]
    filled_relative_mass = ase.data.atomic_masses[filled_atom_number]
    if is_molecule:
        filled_atom_number_2 = ase.data.atomic_numbers[filled_element[1]]
        filled_relative_mass_2 = ase.data.atomic_masses[filled_atom_number_2]
        filled_relative_mass = (filled_relative_mass + filled_relative_mass_2)/2

    
    total_volume = stru.get_volume()

    # calculate the filled volume -> gived volume - occupied volume of existed atoms 
    filled_volume = ((end - start) / cell_length[2]) * total_volume - np.sum((start/cell_length[2] < stru.get_scaled_positions()[:,2]) & (stru.get_scaled_positions()[:,2] < end/cell_length[2])) * np.pi * 0.7 * 0.7
    
    filled_density_A = (filled_density * 6.02214)/(10*filled_relative_mass)  # 个数/埃
    filled_number = np.ceil(filled_density_A*filled_volume)
    
    filled_atoms = Atoms(cell=cell, pbc=True)
    tries = 0
    max_tries = filled_number * 50  # 设置一个尝试的最大次数以避免无限循环

    start_fraction = start/cell_length[2]
    end_fraction = end/cell_length[2]

    while len(filled_atoms) < filled_number:
        point_fraction = generate_random_point(start_fraction, end_fraction)
        point = np.dot(point_fraction, cell)
        filled_atoms.append(Atom(filled_element[0], point))

        is_appro_pos = is_appropriate_pos(filled_atoms, stru, min_distance)
        
        if not is_appro_pos:
            filled_atoms.pop()
        
        if is_appro_pos and is_molecule:
            # random create the another atom which connected to the first.
            
            tries += 1
            inner_tries = 0
            
            while True:
                inner_tries += 1
                
                x_max = bond_length
                x = np.random.uniform(point[0]-x_max, point[0]+x_max)
                y_max = np.sqrt(np.square(bond_length) - np.square(point[0]-x))
                y = np.random.uniform(point[1]-y_max, point[1]+y_max)
                z_ = np.sqrt(np.square(bond_length) - np.square(x-point[0]) - np.square(y-point[1]))
                z = np.random.choice([point[2]-z_, point[2]+z_])
                
                filled_atoms.append(Atom(filled_element[1], (x,y,z)))
                is_appro_pos_2 = is_appropriate_pos(filled_atoms, stru, min_distance, True)
                
                if not is_appro_pos_2:
                    filled_atoms.pop()
                else:
                    break
                
                if inner_tries > 50:
                    filled_atoms.pop()
                    break

        tries += 1
        if tries >= max_tries:
            print("achieve the max tries, please retry or reduce the atoms density or minimal distance between the atoms.")
            for atom in filled_atoms:
                stru.append(atom)
            write("checkPoint.vasp", stru, format="vasp")
            return False
        
        # progress
        n_p = int((len(filled_atoms)/filled_number)*50)
        print("progress: ","*"*n_p, f"{n_p * 2}%; filled: {len(filled_atoms)} (total: {filled_number} {'-'.join(filled_element)}); loop: {tries} (max tries: {max_tries})", end="\r")

    for atom in filled_atoms:
        stru.append(atom)
        
    if prm.visualize==True:
        view(stru)
    
    write("out.vasp", stru, format="vasp")

    print(f"\nEnd of calculation. total {tries} tries.\n")
    
    return True

# Starting
starttime = time.time()
print("Author: gchen")
print("Starting calculation at", end=' ')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))

main()

endtime = time.time()
runtime = endtime-starttime
print(time.strftime("%H:%M:%S on %A %d %B %Y"))
print("Program was running for %.2f seconds." % runtime)