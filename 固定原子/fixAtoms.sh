#!/bin/bash
#使用方法
#0.把脚本和POSCAR放在一个文件夹内
#1.命令行输入"脚本名 固定数目", eg fixAtoms 5, 表示沿着c轴固定最底下5个原子


fix_num=$1 # 固定数目
####################
if [[ ${fix_num} == "" ]];then
	echo "请输入参数"
	exit 1
fi
# 转换POSCAR为linux格式
dos2unix POSCAR
# 将POSCAR中的坐标放入list文件中
sed -nr -e '/^Dir|^Car/,/^\s*$/p' POSCAR | sed -r '/^Dir|^Car|^\s*$/d' | awk '{print NR,$1,$2,$3}' |sort -t' ' -k4 > list
# 将POSCAR中坐标前的文件单独存放
sed -r '/Dir|Car/iSelective Dynamics' POSCAR | awk 'NR==1,/Dir|Car/{print $0}' > POSCAR_HEAD
count=1
# 在遍历POSCAR每一行坐标之前先要改变IFS
IFS_OLD=$IFS
IFS=$'\n'
for row in `cat list`
# 固定原子-->一行一行读取
do
	if ((count <= fix_num)); then
		new_row=$row" F F F"
	else
		new_row=$row" T T T"
	fi
	sed -i ${count}s/${row}/${new_row}/g list
	((count++))
done
# 循环结束将IFS改回来
IFS=$IFS_OLD
#将坐标按照原来顺序排序，并且删除行号，连接POSCAR首尾放入指定文件夹
cat POSCAR_HEAD > POSCAR
sort -n -t' ' -k1 list | cut -d' ' -f2,3,4,5,6,7 >> POSCAR 

rm list
rm POSCAR_HEAD
