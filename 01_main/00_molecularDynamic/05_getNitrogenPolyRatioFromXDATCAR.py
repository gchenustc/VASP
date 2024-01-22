"""
Desc: 
    1. This script calculates the polymerization degree of N atoms for each step in molecular dynamics. It requires the XDATCAR and POSCAR files and returns the polymerization degree at each step, along with a text documentation and a corresponding figure. The polymerization degree is defined as the number of polymerized N atoms divided by the total number of N atoms. An N atom is considered polymerized if it is connected to two other N atoms.
    2. By default, multiprocessing is disabled in the script. If enabled, it only has an effect on the local machine and does not work with the slurm system (such as pbs). Please note that the script has not been tested with other systems.
Author: 
    gchen
Time: 
    2024-1-18 
"""

from ase.io.vasp import read_vasp_xdatcar
from ase.io import read
from ase.geometry.analysis import Analysis
import time
import multiprocessing
import matplotlib.pyplot as plt
import sys

##### setting ######
processing = 1
step = int(sys.argv[1])
##### setting ######

xdatcar = read_vasp_xdatcar('XDATCAR', index=slice(-1))
#poscar = read("POSCAR",format="vasp")
#xdatcar.insert(0,poscar)
xdatcar = xdatcar[:step] #取前10000个构型

atomic_symbols = xdatcar[0].get_chemical_symbols()
total_N = sum(i=="N" for i in atomic_symbols)

def add_empty_list(dictionary, key):
    """
    If the key exists in the dictionary or the value of the key is not empty, return the dictionary as it is.
    If the key does not exist in the dictionary or the value of the key is empty, create the key and an empty list as its value.
    """
    if key not in dictionary or not dictionary[key]:
        dictionary[key] = []
    return dictionary

def get_poly_rato(stru):
    """
    Calculate the polymerization degree of N atoms in a given structure.
    The polymerization degree is defined as the number of polymerized N atoms divided by the total number of N atoms.
    An N atom is considered polymerized if it is connected to two other N atoms.
    """
    bonds_dict = {}  # {1:[2,3],....} # The atoms bonded to atom with index 1 are 2 and 3
    ana = Analysis(stru)
    bonds = ana.get_bonds("N","N")[0] # [(1,2),(1,3),...]
    for idx1,idx2 in bonds:
        add_empty_list(bonds_dict, idx1)
        add_empty_list(bonds_dict, idx2)
        bonds_dict[idx1].append(idx2)
        bonds_dict[idx2].append(idx1)

    count = sum(len(bonds_dict[key]) >= 2 for key in bonds_dict)
    return count/total_N


def write_file(name, lst):
    """
    Write the elements of a list to a file, with each element on a new line.
    """
    with open(name, 'w') as file:
        # Write the content line by line
        for idx,i in enumerate(lst):
            file.write(f"{str(idx+1)}\t")
            file.write(f"{str(i)}\n")


def plot(name, lst, xlabel,ylabel):
    """
    Plot a graph using the given list of values.
    """
    plt.plot(range(1,len(lst)+1), lst)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(name)


def main():
    """
    The main function that executes the script.
    """
    start_time = time.time()  # Record the start time

    with multiprocessing.Pool(processing) as pool:
        ratio_list = list(pool.map(get_poly_rato, list(xdatcar)))

    end_time = time.time()  # Record the end time
    print("Calculation completed in", end_time - start_time, "seconds")

    write_file('polyRatio.txt', ratio_list)
    plot("polyRatio.png", ratio_list, 'step', 'polyDegree')

if __name__ == "__main__":
    main()
