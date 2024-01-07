#!/bin/env bash
# Description: giving the specific parameters to create correspoing name's folders where calculating the phonon, needs POSCAR(should change names, the POSCARs' names are as same as the paremeters) INCAR vasp_job.sh SuperCell.py in current folder. Also needs to install the vaspkit.
# Author: gchen
# Time: 2023-12-27 21.30

supercell="2 2 2"
for i in "$@"
do

    mkdir ${i}_
    cp ${i} ${i}_/POSCAR
    cp vasp_job.sh ${i}_
    cp SuperCell.py ${i}_
    cp INCAR ${i}_

    cd ${i}_

    python SuperCell.py -d --dim ${supercell}
    mv POSCAR POSCAR.bak
    mv POSCAR_super.vasp POSCAR
    echo -e "103" | vaspkit
    sbatch vasp_job.sh

    cd ..

done
