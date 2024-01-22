#/bin/env bash
#desc: move the destinate file to the current folder. the first parameter is the file name and the second parameter is the file depth (referrring to the position of the target file within the current folder hierarchy). 

function judgeLevel(){
    if (( "$1" == 0 ));then
        return 1
    fi
    return 0
}

level=$2

if [ $level -le 0 ];then
    echo "the second parameter need bigger than 0"
    exit 1
fi

mapfile -t arr < <(find . -mindepth $level -maxdepth $level  -name $1 -type f)

mkdir $1
echo "the number of total destinated file: ${#arr[@]}"
echo "the file path are repectively:"
for path in ${arr[@]};
do
    tmp=$(echo ${path:2} | sed 's/\//-/g')
    cp ${path} $1/$tmp
    echo ${path}
done
