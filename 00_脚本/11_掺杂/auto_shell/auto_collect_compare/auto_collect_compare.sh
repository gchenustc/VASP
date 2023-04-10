#!/bin/env bash
# des: 将子目录的CONTCAR全部收集到${work_dir}中，再逐一对比与结构${tem_stru}的相似度（只考虑格点而不考虑原子类型)，输出保存在${log}中，搭配结构收集脚本和结构相似性检测脚本一起

# 用来对比的结构名
tem_stru=BPN_30GPa_3x1x2.vasp
log="out.log"
work_dir="CONTCARs"

# 收集CONTCAR
bash collectContcar.sh ${work_dir}
#for i in `ls`;
#do
#    if [ -f ${i} ] && [[ ${i} =~ ".vasp" ]];then
#        python SuperCell.py -c ${i} -d --dim "2 2 2"
#        rm ${i}
#    fi
#done

cd ${work_dir}
for i in `ls`
do
    if [[ ${i} =~ ".vasp" ]];then
        echo -n "Are >>> ${i} <<< not decompose?    " >> ../${log}
        out=`python ../StruEquiTest.py  ../${tem_stru} ${i}`
        if [ ${out} == "True" ];then
            echo ">>>  True  <<<"  >>  ../${log}
        else
            echo False >> ../${log}
        fi
    fi
done
