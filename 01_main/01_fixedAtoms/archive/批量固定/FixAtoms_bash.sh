#!/bin/bash
############使用方法###############
#遍历当前目录下的子目录，进入子目录下固定POSCAR
##########设定固定数目##########
fix_num=16
####################
for i in `dir`
do
if [ -d ${i} ];then
cd ${i}
# 将POSCAR中的坐标放入list文件中
sed -nr -e '/^Dir|^Car/,/^\s*$/p' POSCAR | sed -r '/^Dir|^Car|^\s*$/d' | awk '{print NR,$1,$2,$3}' |sort -n -t' ' -k4 > list
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
cd ..
fi
done
