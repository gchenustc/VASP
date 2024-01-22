"""
Desc: 
    This script calculates the energy (E0) at each step in molecular dynamics for VASP. It requires the OSZICAR file and returns the energy at each step, along with a text documentation and a corresponding figure. The first param is the number of steps to extract for.
author: 
    gchen
time: 
    2024-1-18 
"""

import matplotlib.pyplot as plt
import subprocess
import sys

def extract_energy(step):
    """
    Extracts energy values from the OSZICAR file and saves them in a text file.
    
    Args:
        step (int): The number of steps to extract energy values for.
    
    Returns:
        None
    """
    # Use awk command to extract energy values from the OSZICAR file and write them to out.txt file
    subprocess.run('awk -F" " \'{if($8=="E0="){printf "%s\\n", $9}}\' OSZICAR > out.txt', shell=True)  
    
    # Read energy values from out.txt file and store them in energyList
    with open('out.txt', 'r') as file:
        energyList = [float(line.strip()) for line in file.readlines()][:step]
    
    # Write energy values to a file
    with open('energy.txt', 'w') as file:
        for idx, energy in enumerate(energyList):
            file.write(f"{idx+1}\t{energy}\n")
    
    # Plot the energy values using matplotlib
    plt.plot(range(1, len(energyList)+1), energyList)
    plt.xlabel('step')
    plt.ylabel('energy')
    plt.savefig("energy.png")
    
    # Delete the out.txt file    
    subprocess.run("rm out.txt", shell=True)

step = int(sys.argv[1])
extract_energy(step)
