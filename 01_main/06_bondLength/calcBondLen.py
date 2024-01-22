"""
Desc: Calculate bond lengths and average bond lengths from CONTCAR file using ase and pymatgen packages.
Usage: Run the script with "python script.py POSCAR"

Requires the following packages:
- ase
- pymatgen

Author: gchen
Time: 2023
Version: 2.0
"""
from ase.io import read
from ase.geometry.analysis import Analysis
from collections import Counter
from itertools import combinations_with_replacement
import sys

# param for structure's name
stru_name = sys.argv[1]

# Read the CONTCAR file using ase
stru_ase = read(stru_name, format="vasp")

# Create an Analysis object from the ase structure
ana = Analysis(stru_ase)

# Get the chemical elements and their counts from the formula
chemSymbolsList = stru_ase.get_chemical_symbols() # ["N", "N", "F", "F", ...]
kinds = Counter(chemSymbolsList) # {'P': '1', 'N': '7'}

# Initialize dictionaries to store bond indices and bond values
bonds_value = {}

# Iterate through each pair of elements and calculate bond values
comb = combinations_with_replacement(kinds.keys(),2)
for group in list(comb): 
    kind1 = group[0]
    kind2 = group[1]
    bond_index = ana.get_bonds(kind1, kind2, unique=True)
    if bond_index[0]:  # Check if there are bonds between kind1 and kind2
        bonds_value[f"{kind1}-{kind2}"] = [bond_index[0], ana.get_values(bond_index)[0]]

# Print the bond values for each pair of elements
for key, value in bonds_value.items():
    print(f"{key}'s bond indices: {value[0]}\n")
    print(f"{key}'s bond lengthes: {', '.join(str(item) for item in value[1])}\n")

print()

# Calculate and print the average bond length for each pair of elements
for key, value in bonds_value.items():
    average = sum(value[1]) / len(value[1])  # Calculate the average value
    print(f"The average of {key} is {average}")