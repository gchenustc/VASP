"""
Desc: 
    This script calculates the temperature at each step in molecular dynamics for VASP. It requires the OSZICAR file and returns the temperature at each step, along with a text documentation and a corresponding figure. The first param is the number of steps to extract for.
author: 
    gchen
time: 
    2024-1-18 
"""

import matplotlib.pyplot as plt
import subprocess
import sys

def extract_temperature(step):
    """
    Calculates the temperature at each step in molecular dynamics for VASP.

    Args:
        step (int): The number of steps to extract for.

    Returns:
        None. Generates a text file 'temperature.txt' with the temperature at each step,
        and a corresponding figure 'temperature.png' showing the temperature vs. step.

    Author: gchen
    Time: 2024-1-18
    """

    # Extract energy values from OSZICAR file
    subprocess.run('awk -F" " \'{if($8=="E0="){sub("\\.$", "", $3);printf "%s\\n", $3}}\' OSZICAR > out.txt', shell=True)  

    file_path = 'out.txt'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        temList = [float(line.strip()) for line in lines][:step]

        # Write temperature values to file
        with open('temperature.txt', 'w') as file:
            for idx, i in enumerate(temList):
                file.write(f"{str(idx+1)}\t")
                file.write(f"{str(i)}\n")

        # Plot temperature vs. step using matplotlib
        plt.plot(range(1, len(temList)+1), temList)
        plt.xlabel('step')
        plt.ylabel('temperature')
        plt.savefig("temperature.png")
        # Remove temporary file
        subprocess.run("rm out.txt", shell=True)

step = int(sys.argv[1])
extract_temperature(step)