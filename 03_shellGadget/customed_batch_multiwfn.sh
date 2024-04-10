# !/bin/env bash
# DESC: Create the script for operating the Multiwfn, the created script need give a parameter(represent the file will operate).
# Time: 20240409
# Author: gchen

echo "please input the script name:"
read -r script_name

echo -e "Please write the multiwfn command, if want to end, input the -> quit \n\$1 refers to the first parameter (file name) for generated script. \n^1 refers to the dirname for first parameter of generated script. \n^2 refers to the basename for first parameter of generated script. \n..."
while IFS= read -r line; do
    if [ "$line" == "quit" ]; then
        break
    fi  
    user_input+=("$line")
done

echo "# !/bin/env bash" > $script_name
echo "# Desc: the Automatic generated script by $0" >> $script_name
echo "# Time: $(date +%F)" >> $script_name
echo 'Multiwfn $1 << EOF >/dev/null' >> $script_name

for i in ${user_input[@]};do

    dirname="\$(echo \$1 | awk -F'.' '{print \$1}')"
    basename="\$(echo \$1 | awk -F'.' '{print \$2}')"
    result=$(echo "$i" | sed "s/\^1/$dirname/g" | sed "s/\^2/$basename/g")                                                                  
    echo $result >> $script_name
done
echo "EOF" >> $script_name

chmod u+x $script_name

echo "Done..."
