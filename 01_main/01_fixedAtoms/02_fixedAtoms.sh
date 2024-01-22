#!/bin/bash
#Desc: Fixing Atoms for the POSCAR format file.
#Usage: ./script.sh 16 # fixing the 16 atoms which the distance from low to high of z axis.
#Author: gchen

fix_num=$1 # Set the number of atoms to fix
if [[ ${fix_num} == "" ]];then
	echo "Please enter the parameter"
	exit 1
fi

# Convert POSCAR to Linux format
dos2unix POSCAR

# Extract coordinates from POSCAR and store them in a list file
sed -nr -e '/^Dir|^Car|^car|^dir/,/^\s*$/p' POSCAR | sed -r '/^Dir|^Car|^car|^dir|^\s*$/d' | awk '{print NR,$1,$2,$3}' | sort -n -t' ' -k4 > list

# Separate the coordinates from the original POSCAR file
sed -r '/Dir|Car|dir|car/iSelective Dynamics' POSCAR | awk 'NR==1,/Dir|Car|dir|car/{print $0}' > POSCAR_HEAD

count=1

# Iterate through each row of the list file
# Fix atoms one by one
IFS_OLD=$IFS
IFS=$'\n'
for row in `cat list`
do
	if ((count <= fix_num)); then
		new_row=$row" F F F"
	else
		new_row=$row" T T T"
	fi
	sed -i ${count}s/${row}/${new_row}/g list
	((count++))
done

# Restore the original IFS value
IFS=$IFS_OLD

# Sort the coordinates in the original order and remove line numbers
# Append the sorted coordinates to the end of the POSCAR file
cat POSCAR_HEAD > POSCAR
sort -n -t' ' -k1 list | cut -d' ' -f2,3,4,5,6,7 >> POSCAR 

# Remove temporary files
rm list
rm POSCAR_HEAD
