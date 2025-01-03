#!/bin/env sh
# Author: gchen
# Description: This script processes and compares the results of DFT and DP (Deep Potential) models before and after training for multiple systems. It creates directories for each iteration and system, runs the DP test, and generates plots using a Python script.
# Parameter: 
#   $1: The number of systems to process
#   $2: The current iteration number (integer)
# Usage: 
#   ./script_name.sh <number_of_systems> <current_iteration>

f_py_plot="plot_error_DP_DFT.py"

# 检查参数是否存在
if [ $# -ne 2 ]; then
  echo "Usage: $0 <number of systems> <iteration: int>"
  exit 1
fi

# 检查绘图文件是否存在
if [ ! -e "$f_py_plot" ]; then
  echo "Need plotting py file: ${f_py_plot}"
  exit 1
fi

n_sys=$1
cur_iter=$2

index_sys=$(( n_sys - 1 ))
after_iter=$(( cur_iter + 1 ))

mkdir iter.${cur_iter}

cd iter.${cur_iter}

for i in $(seq 0 $index_sys)
do

    echo -e "/n****** current is the data.${i} *****/n"

    if ! test -d "../../iter.00000${cur_iter}/02.fp/data.00${i}";then
        echo -e "/n****** data.${i} is not existed *****/n"
        continue
    fi

    mkdir sys_${i}
    cd sys_${i}

    mkdir pre_train
    mkdir after_train

    cd pre_train
    cp ../../../${f_py_plot} ./

    echo -e "/n****** current is pre_train results of the data.${i} *****/n"

    dp test -m ../../../../iter.*${cur_iter}/00.train/000/frozen_model.pb -s ../../../../iter.*${cur_iter}/02.fp/data.00${i}  -l out.log   -v 2 -d ./detile
    python ${f_py_plot}
    cd ..

    cd after_train
    cp ../../../${f_py_plot} ./

    echo -e "/n****** current is after_train results of the data.${i} *****/n"

    dp test -m ../../../../iter.*${after_iter}/00.train/000/frozen_model.pb -s ../../../../iter.*${cur_iter}/02.fp/data.00${i}  -l out.log   -v 2 -d ./detile
    python ${f_py_plot}
    cd ..

    cd ..

done

cd ..

