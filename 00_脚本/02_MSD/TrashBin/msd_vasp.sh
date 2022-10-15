#!/bin/env bash
#Describe 用于在vasp中计算msd，需要POSCAR和XDATCAR文件
#Time 2022-1-11
#Author gchen


#提取POSCAR中的晶格常数>cons.csv
sed -n '3,5p' POSCAR | tr -s ' ' | sed -nr 's/ (.*)/\1/p' > cons.csv
#这串代码也可以 sed -n '3,5p' POSCAR | tr -s ' ' | sed -r 's/.//' > cons.csv

#提取POSCAR中的初始位置>init_POS.csv
sed -n '/Direct/,$p' POSCAR| grep -v 'Direct' | tr -s ' ' | sed 's/.//' > init_pos.csv

#提取原子信息>atoms_info.csv
sed -n '6,7p' POSCAR |tr -s ' '| sed 's/.//' > atoms_info.csv

#提取XDATCAR的位置信息>xdatcar
awk 'NR>7' XDATCAR | grep -v 'Direct' | tr -s ' ' |sed 's/.//' > xdatcar.csv

#获得步数和步长>xdatcar
grep 'Direct' XDATCAR|wc -l > step.csv
read -p "输入步长" len
echo $len >> step.csv

python msd_vasp.py

rm cons.csv
rm init_pos.csv
rm atoms_info.csv
rm xdatcar.csv
rm step.csv


