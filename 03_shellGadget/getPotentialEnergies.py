import subprocess
import os
import sys
from ase.io.vasp import read_vasp_out
from ase.io import read
import numpy as np
"""
Desc: This script retrieves potential energy in specified subfolders(appointed by params) and appends the data to a file named potential_energy.txt.
Author: gchen
Time: 20240122
"""

# folder lists; the process will enter these folder and execute.
folder_list = sys.argv[1:]

# write the title for log
with open('potential_energy.log', 'w') as log_file:
    log_file.write(subprocess.check_output(['date']).decode('utf-8').strip() + '\n')
    log_file.write("folder\tpotential_energy\tpotential_energy/atom\n")

info_list = []
# enter the appointed folders and obtain forces.
for folder in folder_list:
    if os.path.isdir(folder):
        #os.chdir(folder)
        atoms = read(f'{folder}/OUTCAR', format='vasp-out')
        outcar= read_vasp_out(file=f'{folder}/OUTCAR', index=-1)
        energy = outcar.get_potential_energy()
        average_energy = energy/len(atoms)
        info_list.append(f"{folder}\t{energy}\t{average_energy}\n")
        #os.chdir("..")

# write the forces information
with open('potential_energy.log', 'a') as log_file:
    log_file.writelines(info_list)
