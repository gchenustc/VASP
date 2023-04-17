#!/bin/env bash
# 将子目录下的文件夹全部移到(不包含名为attach的文件夹) ../../下

for i in {0..4..1}
do
    cd ${i}

    for i in `ls`;
    do
        if [ ${i} != attach ] $$ [ -d ${i} ];then
            mv ${i} ../../
        fi
    done

    cd ..
done
