#!/bin/bash
#Desc: In this script, we have a batch processing shell script that adds Python files from a specified folder to the .bashrc file as aliases. Additionally, it checks whether the alias for a specific file already exists in the .bashrc file and skips adding it if it does. The script ensures that the .bashrc file is reloaded to immediately apply the changes.
#Author: gchen
#Time: 20240119

# Specify the folder path
folder="/public/home/c1337375425/scripts"

# Iterate through the Python files in the folder
for file in $folder/*.py; do
    # Extract the filename (without the path and extension)
    filename=$(basename "$file" .py)

    # Check if the alias already exists in the .bashrc file
    if ! grep -q "alias ${filename}=" ~/.bashrc; then
        # Add the alias to the .bashrc file
        echo "alias ${filename}=\"python ${file}\"" >> ~/.bashrc
        echo "Added alias for ${filename}"
    else
        echo "Alias for ${filename} already exists, skipping"
    fi
done

# Reload the .bashrc file
source ~/.bashrc

echo "Aliases have been added to the .bashrc file"

