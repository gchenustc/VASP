#!/bin/env bash
# Desc: auto submit the vasp task, the current folder needs the POSCAR file. Noted, you should to assagin the INCAR and submitted script path and edit them in the script
# Usage: ./xx.sh


echo "please select a INCAR file: "
echo `ls ~/storage/vasp/`
read choice

cp ~/storage/vasp/$choice ./INCAR
cp ~/scripts/sub_scripts/vasp_job.sh ./

vaspkit << EOF >> /dev/null
103
EOF

echo "Do you want to fix atoms? (y/n)"
read choice

if [ "$choice" = "y" ]; then

    echo "You chose to proceed."
    cp ~/scripts/fixAtoms.sh ./

    echo "input the fixed number"
    read number
    ./fixAtoms.sh $number

elif [ "$choice" = "n" ]; then
    echo "You chose not to proceed."

else
    echo "Invalid choice. the fix operation not proceed."
fi

echo "Do you want to submit task? (y/n)"
read choice

if [ "$choice" = "y" ]; then

    echo "You chose to proceed."
    sbatch vasp_job.sh

elif [ "$choice" = "n" ]; then
    echo "You chose not to proceed."

else
    echo "Invalid choice. the submit operation not proceed."
fi
