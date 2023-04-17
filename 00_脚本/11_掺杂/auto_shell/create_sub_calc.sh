#!/bin/env bash
# 生成文件夹名从1到5的文件夹，将attach文件夹和batch_calc.sh拷贝进每一个，再依次运行每个文件夹下的batch_calc.sh

for i in {1..5..1}
do
    cp -r attach ${i}
    cp batch_calc.sh ${i}

    cd ${i}
    bash batch_calc.sh
    cd ..

    sleep 10000
done
