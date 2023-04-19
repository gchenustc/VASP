#!/bin/env bash
# desc: print sub-folder's enthalpy and compile them into the enthalpy.txt文件夹
# author: gchen
echo "time: `date`" > enthalpy.txt

for i in `ls` ;do
    if [ -d ${i} ];then
        cd ${i}
        echo $(awk -F" " -v label=${i} 'BEGIN{OFS="\t"};{if($4=="E0="){print label,$5}}' OSZICAR | tail -1) >> ../enthalpy.txt
        cd ..
    fi
done
