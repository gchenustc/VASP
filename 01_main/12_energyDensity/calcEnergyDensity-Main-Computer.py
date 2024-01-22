"""
Desc:
    This script calculates the energy density of a given chemical reaction based on the ground state combination,ratios, and energies of the involved elements.
    For example, for the reaction 2CO + O2 = 2CO2,
    set ground_state_combinationination = [{"C": 1, "O": 1}, {"O": 2}]
    ratio = [2, 1]
    energy = [1 * energy_of_CO, 1 * energy_of_O2] : Corresponding energies for the ground state combination
Author:
    gchen

Energy values for the ground state elements(eV)
P_Fibrous_0GPa = -5.54817
N_alpha_0GPa = -8.37770
N_epsilon_0GPa = -8.37539
Si_0GPa = -5.424965
"""

from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.io.vasp import Poscar
from pymatgen.core.periodic_table import Element
from collections import Counter
from ase.io import read


ground_state_combination = [{"N": 1}, {"Si": 1}]  # Number of each element in the ground state combination
ratio = [6, 1]  # Ratios of the elements in the ground state combination
energy = [-8.37539, -5.424965]  # Energies of the ground state elements
# Editable parts

# Multiply the energies by the corresponding ratios
for i in range(len(energy)):
    energy[i] = energy[i] * ratio[i]

# Read energy from OSZICAR
oszicar = Oszicar(filename="OSZICAR")
final_energy_product = oszicar.final_energy

# Read structure from POSCAR
poscar = Poscar.from_file(filename="POSCAR")
stru = poscar.structure
stru_ase = read("POSCAR", format="vasp")

# Get chemical symbols
symbols = stru_ase.get_chemical_symbols()  # ["F", "N", "N", ...]

# Count the occurrences of each symbol
symbols_sta = Counter(symbols)  # Counter({'F': 8, 'N': 8})

ground_state_sta = Counter()  # Counter for the ground state {'N': 3, 'F': 3}
for i, item in enumerate(ground_state_combination):
    for key, value in item.items():
        item[key] = value * ratio[i]
    ground_state_sta.update(item)

# Check if the input is correct
keys = sorted(list(set(symbols_sta.keys()).union(set(ground_state_sta.keys()))))  # ["F, "N"]
ratio_ = symbols_sta[keys[0]] / ground_state_sta[keys[0]]
for key in keys:
    if symbols_sta[key] / ground_state_sta[key] != ratio_:
        print("Incorrect ground state ratio")
        raise KeyError

# Determine the quantity of the final state substance
n = symbols_sta[keys[0]] / ground_state_sta[keys[0]]

# Calculate energy density
# Numerator
numerator = final_energy_product / n
for i in range(len(ratio)):
    numerator -= energy[i]

mass = sum(
    ground_state_sta[key] * Element(key).atomic_mass.real for key in keys
)
print("Energy density:", numerator * 96.485 / mass)