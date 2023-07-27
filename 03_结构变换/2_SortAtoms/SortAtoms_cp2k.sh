#!/bin/env bash
# 作者：陈果
# 描述：cp2k原子排序脚本,按照z轴从小到大排列


# ####### 输出文件名 #########
input=subsys.inc

# 将原子位置拷贝到tmpfile
grep -P "[\s]*N([+\s]+[\d]+\.[\d]*){3}" ${input} > tmpfile

# 将tmpfile的z轴进行排序，并拷贝到sorted_tmpfile中
awk '{print $NF,$0}' tmpfile | sort -g | cut -d" " -f2- > sorted_tmpfile  # -f2- 排列z轴

# 将sorted_tmpfile中的内容替换输入文件的原子位置
sorted_line=1
while read line; do
  # 获得当前${line}所在行的行号 && 将找到的${line}的上一行插入sorted_tmpfile文件的${sorted_line}行，该行依次递增 && 删除原来的行
  n_line=$(sed -n "/${line}$/=" ${input})
  sed -i "/${line}$/i \  $(sed -n "${sorted_line} p" sorted_tmpfile)1" ${input}
  (( n_line+=1 ))
  sed -i "${n_line} d" ${input}
  (( sorted_line+=1 ))
done < tmpfile


rm tmpfile
rm sorted_tmpfile

