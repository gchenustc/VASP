#/bin/env bash
#desc: check if achieve convergences of vasp relax calculations for sub-folders. the first param is the folder depth and the second param is optional, if exist, the programme only check the folder has the name of second parma.

if [[ "$1" -le 0 ]];then
    echo "the first parameter need bigger than 0"
    exit 1
fi

if [ $# -eq 2 ]; then
    mapfile -t arr < <(find . -mindepth $1 -maxdepth $1 -name $2 -type d)
    echo 1
else
    mapfile -t arr < <(find . -mindepth $1 -maxdepth $1 -type d)
fi

echo "the number of total folder: ${#arr[@]}"
echo "the relax results  are repectively:"
for path in ${arr[@]};
do
    echo -n "${path}: "
    cd ${path}
    echo `grep reached OUTCAR|tail -1`
    cd - > /dev/null
done
