import subprocess
import os
import sys
from ase.io.vasp import read_vasp_out
import numpy as np
"""
Desc: appointed the params which are the folder' names for the Vasp package, the proceedings will enter these folders and get each calculation's forces and store them in the current folder's file.
Author: gchen
Time: 20240119
"""

# folder lists; the process will enter these folder and execute.
folder_list = sys.argv[1:]

# write the title for log
with open('forces.log', 'w') as log_file:
    log_file.write(subprocess.check_output(['date']).decode('utf-8').strip() + '\n')
    log_file.write("folder\tmax_force\tave_force\n")

info_list = []
# enter the appointed folders and obtain forces.
for folder in folder_list:
    if os.path.isdir(folder):
        os.chdir(folder)

        outcar= read_vasp_out(file='OUTCAR', index=-1)
        forces = outcar.get_forces()
        _norm = np.linalg.norm(forces,axis=1)
        max_forces = _norm.max()
        ave_forces = _norm.mean()
        
        info_list.append(f"{folder}\t{max_forces:.6f}\t{ave_forces:.6f}\n")

        os.chdir("..")

# write the forces information
with open('forces.log', 'a') as log_file:
    log_file.writelines(info_list)
