#!/bin/bash
#DESC：This script converts a geometry file (e.g., .xyz, .pdb, .inp etc.) into different output formats using Multiwfn. The supported output formats are POSCAR, CIF, PDB, and XYZ. The script prompts the user to choose the desired output format and then uses Multiwfn to perform the conversion.
# USAGE: ./script.sh input_file.inp/vasp/cif

if [ "$1" ]; then
    input_file="$1"
else
    echo -e "\033[31mERROR:\033[0m Missing file operand! Please identify the name of the geometry file!"
    exit 1
fi

echo -e "Choose the format of the output file:\n1. POSCAR \n2. cif \n3. pdb\n4. xyz"
read -p "Enter your choice: " output_format_choice

case $output_format_choice in
    1)
        output_format=27
        ;;
    2)
        output_format=33
        ;;
    3)
        output_format=1
        ;;
    4)
        output_format=2
        ;;
    *)
        echo -e "\033[31mERROR:\033[0m Invalid choice!"
        exit 1
        ;;
esac

echo -e "$input_file\n100\n2\n$output_format\n\n0\nq" | Multiwfn 1>>/dev/null

echo "Done!"
