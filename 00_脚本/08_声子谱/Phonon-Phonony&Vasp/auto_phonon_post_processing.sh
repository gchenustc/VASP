#!/bin/emv bash

echo -e "305\n3"| vaspkit
mv KPATH.phonopy band.conf
sed -i "s/DIM =/DIM = 3 1 3/1" band.conf  # !!!! attention it's need to change dim
sed -i "s/FORCE_CONSTANTS = READ/FORCE_CONSTANTS = WRITE/1" band.conf
phonopy -f ./disp-*/vasprun.xml
phonopy -c POSCAR band.conf -p -s
phonopy-bandplot --gnuplot > phonon.out
