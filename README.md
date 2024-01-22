# VASP SCRIPTS

Please refer to the script for detailed usage instructions. The following content is only an introduction to the script.

## 01_main
### 01_main/00_molecularDynamic
#### 01_main/00_molecularDynamic/01_msd/msd.py
This script calculates the mean square displacement (MSD) of atoms in a molecular dynamics simulation and generates plots of the MSD. The MSD measures the average displacement of atoms from their initial positions over time.

#### 01_main/00_molecularDynamic/02_getEnergyFromOSZICAR.py
This script calculates the energy (E0) at each step in molecular dynamics for VASP. It requires the OSZICAR file and returns the energy at each step, along with a text documentation and a corresponding figure. The first param is the number of steps to extract for.

#### 01_main/00_molecularDynamic/03_getTemperatureFromOSZICAR.py
This script calculates the temperature at each step in molecular dynamics for VASP. It requires the OSZICAR file and returns the temperature at each step, along with a text documentation and a corresponding figure. The first param is the number of steps to extract for.

#### 01_main/00_molecularDynamic/04_getBondLenFromXDATCAR.py
This script calculates the appointed bond length for each step in molecular dynamics. It requires the XDATCAR and POSCAR files and returns the appointed bond length at each step, along with a text documentation and a corresponding figure. 

#### 01_main/00_molecularDynamic/05_getNitrogenPolyRatioFromXDATCAR.py
This script calculates the polymerization degree of N atoms for each step in molecular dynamics. It requires the XDATCAR and POSCAR files and returns the polymerization degree at each step, along with a text documentation and a corresponding figure. The polymerization degree is defined as the number of polymerized N atoms divided by the total number of N atoms. An N atom is considered polymerized if it is connected to two other N atoms.

### 01_main/01_fixedAtoms
#### 01_main/01_fixedAtoms/01_fixedAtoms.sh(02_fixedAtoms.py)
Fixing Atoms for the POSCAR format file.

#### 01_main/01_fixedAtoms/03_fixedAtoms-center.sh
Fixing Atoms which at the position of the center for z axis.

### 01_main/02_RDF
#### 01_main/02_RDF/02_rdf.py
Calculating the rdf for molecular dynamics (vasp). the first and the second params are the frame indices for the calculated range. the third param is the max distance for the rdf (rmax). this value less than the half of the shortest axis' length. the fourth param is the number of division for the rmax. the fifth param is the appointed element. 

### 01_main/05_COHP
#### 01_main/05_COHP/01_calcCohp-ase.py
This script allows you to calculate the COHP (Crystal Orbital Hamilton Population) using VASP (Vienna Ab initio Simulation Package) and analyze the results using the Lobster tool. COHP provides insights into the bonding and interactions between atoms in a crystal structure.

#### 01_main/05_COHP/02_plotCohp-ase.py
After calculate COHP, use it.

### 01_main/06_bondLength
#### 01_main/06_bondLength/calcBondLen.py
This script calculates bond lengths and average bond lengths from CONTCAR file using ase and pymatgen packages.

### 01_main/07_band
Calculating band and band's dos using ase (call the VASP calculation). 

### 01_main/08_phonon
Calculating phonon and postprocessing it.

### 01_main/09_structureMatcher
#### 01_main/09_structureMatcher/01_struEquiTest.py
This script detects whether multiple structures are equivalent.

#### 01_main/09_structureMatcher/02_posEquiText.py
This script detect atomic site equivalence.

### 01_main/11_dope
#### ### 01_main/11_dope/dope.py
The provided code is a program that generates doped structures based on a template structure. It replaces specified elements in the template structure with a doping element to create new structures. The program ensures that the generated structures are not duplicated with previously generated structures or existing structures in a given directory.