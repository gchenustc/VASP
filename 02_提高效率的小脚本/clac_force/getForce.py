from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp import Poscar
import argparse
from pymatgen.core.periodic_table import Element  
import ase
from ase.io.vasp import read_vasp_out
import numpy as np
# 计算计算离子步最后一步的最大force和平均force
# 需要OUTCAR

outcar= read_vasp_out(file='OUTCAR', index=-1)
forces = outcar.get_forces()
_norm = np.linalg.norm(forces,axis=1)
max_forces = _norm.max()
ave_forces = _norm.mean()
print(max_forces, ave_forces)

