"""
Author: gchen
Date: 2023.9.8
Desc: This script detect atomic site equivalence. Run this file directly from the command line.
Modify the POSCAR filename and specify the element which to be checked.
"""
from ase.io import write, read
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
import os
import itertools as it

# --- Modify ---
poscar = read("POSCAR", format="vasp")  ## Specify the POSCAR filename
select_kind = "all"  ## Specify the element to compare; "all" checks all atoms; "N" only checks N atoms
# --- Modify ---

def _ase2pymatgen(ase_stru):
    """Converts an ASE structure to a pymatgen structure."""
    write(filename=".vasp", images=ase_stru, format="vasp")
    pymatgen_stru = Structure.from_file(filename=".vasp")
    os.remove(".vasp")
    return pymatgen_stru

def symmetryCheck(atoms1, atoms2):
    atoms1 = _ase2pymatgen(atoms1)
    atoms2 = _ase2pymatgen(atoms2)
    
    # primitive_cell=True: allows comparison of structures with different numbers of atoms; scale: whether to scale to the same volume, disable for high precision
    comp = StructureMatcher(ltol=0.2, stol=0.2, angle_tol=2, primitive_cell=True, scale=False)
    return comp.fit(atoms1, atoms2)

select_idx = []
for idx, atom in enumerate(poscar):
    if select_kind.strip() == "all":
        select_idx.append(idx)
    elif atom.symbol == select_kind:
        select_idx.append(idx)

combi = it.combinations(select_idx, 2)

for idx1, idx2 in combi:
    poscar1 = poscar.copy()
    poscar2 = poscar.copy()
    poscar1[idx1].symbol = "Zr"
    poscar2[idx2].symbol = "Zr"
    print(f"no.{idx1+1}: {poscar[idx1].symbol} ({' '.join(map(str, poscar1[idx1].position))})\
            \n-\n\
no.{idx2+1}: {poscar[idx2].symbol} ({' '.join(map(str, poscar2[idx2].position))})\
            \n-\n\
Equivalence: {symmetryCheck(poscar1, poscar2)}")
    print("-"*15)