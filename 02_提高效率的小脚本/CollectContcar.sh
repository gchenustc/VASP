# 获得子目录中的所有 CONTCAR 并且改名为 dirname_CONTCAR

file=`ls`
for i in $file
do
    if [ -d ${i} ];then
        cd ${i}
		if [ -e CONTCAR ];then
				cp CONTCAR ../${i}_CONTCAR
		fi
		cd ..
    fi
done
        
