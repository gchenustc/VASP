"""
Description: This script detects whether multiple structures are equivalent.
Usage: python *.py POSCAR1 POSCAR2 POSCAR3 # POSCAR* are the paths to the structure files
Author: gchen
"""

from ase.io import write, read
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
import os
import argparse
from itertools import combinations


def _ase2pymatgen(ase_stru):
    """
    Converts an ASE structure to a pymatgen structure.
    
    Args:
        ase_stru: ASE structure object.
    
    Returns:
        pymatgen_stru: Pymatgen structure object.
    """
    write(filename=".vasp", images=ase_stru, format="vasp")
    pymatgen_stru = Structure.from_file(filename=".vasp")
    os.remove(".vasp")
    return pymatgen_stru

def symmetryCheck(atoms1, atoms2):
    """
    Checks the symmetry equivalence of two structures.
    
    Args:
        atoms1: Path to the first structure file.
        atoms2: Path to the second structure file.
    
    Returns:
        results: Boolean value indicating whether the structures are symmetrically equivalent.
    """
    atoms1 = _ase2pymatgen(read(atoms1))
    atoms2 = _ase2pymatgen(read(atoms2))
    
    # primitive_cell=True: allows comparison of structures with different numbers of atoms
    # scale: whether to scale the structures to the same volume (disable for high precision)
    comp = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=1, primitive_cell=True, scale=False)
    return comp.fit(atoms1, atoms2)


def main(stru_list):# keys: combinations of structures to be checked, values: whether they are the same structure (bool)
    for cb in list(combinations(stru_list,2)):
        results = symmetryCheck(cb[0], cb[1])
        print("Returns True if symmetry equivalence exists: ")
        print(f"[{cb[0]}] and [{cb[1]}]: {results}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="None")
    parser.add_argument('POSCAR', nargs="*", help="Specify the POSCAR file")
    prm = parser.parse_args()
    main(prm.POSCAR)