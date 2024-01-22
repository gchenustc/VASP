#!/bin/bash
#Desc: This script retrieves E0 in specified subfolders(appointed by params) and appends the data to a file named enthalpy.txt.
#Author: gchen
#Time: 20240119

# Iterate through the specified folders and perform operations
echo "time: $(date)" > enthalpy.txt
for folder in "$@"; do
    if [ -d "${folder}" ]; then
        cd "${folder}"
        echo $(awk -F" " -v label="${folder}" 'BEGIN{OFS="\t"};{if($4=="E0="){print label,$5}}' OSZICAR | tail -1) >> ../enthalpy.txt
        cd ..
    fi
done

