"""
Desc: 
    1. This script calculates the appointed bond length for each step in molecular dynamics. It requires the XDATCAR and POSCAR files and returns the appointed bond length at each step, along with a text documentation and a corresponding figure. 
    2. By default, multiprocessing is disabled in the script. If enabled, it only has an effect on the local machine and does not work with the slurm system (such as pbs). Please note that the script has not been tested with other systems.
    3. The first parameter is the selected steps. For example: 2000. The second and third parameters are the atom types, and the script will calculate the connected bond length. For example: python script.py 2000 N N.
author: 
    gchen
time: 
    2024-1-18 
"""

from ase.io.vasp import read_vasp_xdatcar
from ase.io import read
from ase.geometry.analysis import Analysis
from collections import Counter
import multiprocessing
import matplotlib.pyplot as plt
import time
import sys

##### setting ######
processing = 1
step = int(sys.argv[1])
bondType = [sys.argv[2],sys.argv[3]] # the bond type stated.
##### setting ######

# Read XDATCAR file
xdatcar = read_vasp_xdatcar('XDATCAR', index=slice(-1))

# Read POSCAR file
poscar = read("POSCAR",format="vasp")

# Insert POSCAR at the beginning of XDATCAR
xdatcar.insert(0,poscar)

# Select the first 'step' configurations
xdatcar = xdatcar[:step]

# Get the chemical symbols of the atoms in the system
chemSymbolsList = xdatcar[0].get_chemical_symbols()

# Count the number of occurrences of each chemical symbol
numCounts = Counter(chemSymbolsList)

# Get the unique chemical symbols
elements = list(set(chemSymbolsList))

# Get the number of unique chemical symbols
numElements = len(numCounts)

def process_data(stru):
    """Calculate the average bond length between two specified atom types"""
    ana = Analysis(stru)
    bondIndexList = ana.get_bonds(bondType[0], bondType[1], unique=True)
    bondLengthList = ana.get_values(bondIndexList)[0]
    return sum(bondLengthList) / len(bondLengthList)

def write_file(name, lst):
    """Write the calculated bond lengths to a file"""
    with open(name, 'w') as file:
        for idx, i in enumerate(lst):
            file.write(f"{str(idx+1)}\t")
            file.write(f"{str(i)}\n")

def plot(name, lst, xlabel, ylabel):
    """Plot the calculated bond lengths"""
    plt.plot(range(1, len(lst)+1), lst)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(name)
    
if __name__ == "__main__":
    start_time = time.time()  # Record the start time

    with multiprocessing.Pool(processing) as pool:
        averBlList = list(pool.map(process_data, list(xdatcar))) 

    end_time = time.time()  # Record the end time
    print("Calculation completed in", end_time - start_time, "seconds")

    write_file('averageNNBondLen.txt', averBlList)
    plot("averageNNBondLen.png", averBlList, 'step', 'average N-N bond length')
