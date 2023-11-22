# 获得子目录中的所有的指定文件（该脚本传入的第一个参数是文件名），移动到当前文件夹中

for i in `ls`
do
    if [ -d ${i} ];then
        cd ${i}
		if [ -e "$1" ];then
				cp "$1" ../${i}_"$1"
		fi
		cd ..
    fi
done
