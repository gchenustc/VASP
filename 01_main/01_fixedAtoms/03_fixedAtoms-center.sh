#!/bin/bash
#Desc: Fixing Atoms for the POSCAR format file. (fixing center atoms)
#Usage: ./script.sh 16 40 # Fix 16 atoms out of 40 atoms in total.
#Author: gchen

fix_num=$1 # Set the number of atoms to fix
atom_num=$2
if [[ ${fix_num} == "" ]];then
	echo "Please enter the parameters"
	exit 1
fi
# Convert POSCAR to Linux format
dos2unix POSCAR
# Extract the coordinates from POSCAR and store them in a list file
sed -nr -e '/^Dir|^Car/,/^\s*$/p' POSCAR | sed -r '/^Dir|^Car|^\s*$/d' | awk '{print NR,$1,$2,$3}' | sort -n -t' ' -k4 > list
# Store the lines before the coordinates in POSCAR separately
sed -r '/Dir|Car/iSelective Dynamics' POSCAR | awk 'NR==1,/Dir|Car/{print $0}' > POSCAR_HEAD

count=1
fix=0
start=$(( (atom_num-fix_num)/2 ))
# Change IFS before iterating through each line of POSCAR
IFS_OLD=$IFS
IFS=$'\n'
# Iterate through each line of the list file
for row in `cat list`
# Fix atoms --> Read line by line
do
	if (( count > start )) && (( fix <  fix_num)); then
		new_row=$row" F F F" # Fix the atom by adding "F F F" at the end of the row
		((fix++))
	else
		new_row=$row" T T T" # Keep the atom movable by adding "T T T" at the end of the row
	fi
	sed -i ${count}s/${row}/${new_row}/g list # Replace the row in the list file with the new row
	((count++))
done

# Restore IFS after the loop ends
IFS=$IFS_OLD
# Sort the coordinates in the original order and remove line numbers, then append them to the end of POSCAR
cat POSCAR_HEAD > POSCAR
sort -n -t' ' -k1 list | cut -d' ' -f2,3,4,5,6,7 >> POSCAR 

rm list
rm POSCAR_HEAD
