# !/bin/env bash
# Desc: convert the .inp or.restart into the .vasp file. Noted, you should install the Multiwfn and the corresponding command is "Multiwfn". The system will auto recongnize the cp2k structure file and convert it.
# Time: 20240409
# Author: gchen

# recongnize the .inp or .restart file and store the name into the array
file_name=($(ls *.inp *.restart))

# define the convert function
convert_file() {
    # get the parameter's content before the .
    dirname=$(echo $1 | awk -F'.' '{print $1}')
    basename=$(echo $1 | awk -F'.' '{print $2}')
    Multiwfn $1 << EOF >/dev/null
    100
    2
    27
    ${dirname}_${basename}.vasp
    0
    q
EOF
}

for i in ${file_name[@]}; do
    echo -e "convert $i ... \n ==============="
    echo
    convert_file $i
done

echo "convert done"
