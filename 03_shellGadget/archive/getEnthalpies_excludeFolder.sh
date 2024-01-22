#!/bin/bash
#Desc: This script get E0 in all subfolders, It also adds the ability to exclude specific folders in appointed params.
#Author: gchen
#Time: 20240119

# Exclude folders specified as arguments
exclude_folders=("$@")

# Function to check if a folder should be excluded
should_exclude() {
    local folder_name="$1"
    for excluded_folder in "${exclude_folders[@]}"; do
        if [ "$folder_name" = "$excluded_folder" ]; then
            return 0  # Folder should be excluded
        fi
    done
    return 1  # Folder should not be excluded
}

# Iterate through the subfolders of the main folder
echo "time: $(date)" > enthalpy.txt
for folder in `ls`; do
    if [ -d ${folder} ];then
        folder_name=$(basename "$folder")
        if should_exclude "$folder_name"; then
            continue  # Skip the folder
        fi

        cd ${folder_name}
        # Enter the folder and perform operations
        echo $(awk -F" " -v label="${i}" 'BEGIN{OFS="\t"};{if($4=="E0="){print label,$5}}' OSZICAR | tail -1) >> ../enthalpy.txt
        cd ..
    fi
done

