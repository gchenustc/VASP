# 获得子目录中的所有 CONTCAR 并且改名为 dirname_CONTCAR.vasp

work_dir=${1} # 收集的CONTCAR所在的目录

if [ -d ${work_dir} ];then
    rm ${work_dir}
    mkdir ${work_dir}
else
    mkdir ${work_dir}
fi


file=`ls`
for i in $file
do
    if [ -d ${i} ];then
        cd ${i}
		if [ -e CONTCAR ];then
				cp CONTCAR ../${work_dir}/${i}_CONTCAR.vasp
		fi
		cd ..
    fi
done
        
