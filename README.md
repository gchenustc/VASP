# VASP脚本手册  
自用的VASP前处理和后处理脚本  
转载时请注明出处  
---
  
# 使用教程
## 功能性脚本
1. Fixatoms.sh  
- 描述：按照沿着z轴的顺序从底部开始固定胞原子  
- 环境：Shell 
- 用法：  
``` shell
chmod u+x Fixatoms.sh
./Fixatoms.sh 16  # 固定POSCAR底部的16个原子
```
---
  
2. Fixatoms.py  
- 描述：按照沿着z轴的顺序从底部开始固定胞原子  
- 环境： python, 需要安装numpy  
``` shell
pip3 install numpy
```
- 用法：  
``` shell
python -c POSCAR -o POSCAR_out -n 10 # 固定POSCAR底部10个原子，输出为POSCAR_out
```
- 批量固定  
``` shell
python -c POSCAR1 POSCAR2 -o POSCAR_out1 POSCAR_out2 -n 10 12 
```
---
  
3. MSDVasp.sh, MSDVasp.py  
- 描述：计算分子动力学的MSD, 只支持立方，四方，正交晶系，因为是早期写的脚本，有些结构可能会出错  
- 环境: shell, python  
python需安装numpy，matplotlib和pandas，如果需要保存图片，要在Linux中创建文件 ~/.config/matplotlib/matplotlibrc。添加 backend : Agg  
``` sheLL
pip install numpy
pip install matplotlib
pip install pands
```
- 用法：  
``` shell
# 1. 将计算完的POSCAR和XDATCAR文件与MSDVasp.sh, MSDVasp.py文件放在同一个目录下
# 2.
chmod u+x MSDVasp.sh
# 3.
./MSDVasp.sh
```
---
  
4. chgdiff.py  
- 描述：计算差分电荷
- 环境: python, 需要安装numpy,ase  
``` shell
pip3 install numpy
pip3 install ase 
```
- 用法：  
``` shell
python chgdiff.py CHGCAR1 CHGCAR2  #输出文件名为 CHGDIFF.
```
- 注意，该脚本来源于SS
> https://github.com/Chengcheng-Xiao/Tools/tree/master/VASP  
使用该脚本请关注 "https://github.com/Chengcheng-Xiao/"  
---
  
5. kpGen.py  
- 描述：产生KPOINTs文件  
- 环境：python, 需要安装numpy,ase, seekpath  
``` shell
pip3 install numpy
pip3 install ase 
pip3 install seekpath 
```
- 用法：  
``` shell
python kpGen.py -c POSCAR  # the most concise usage
python kpGen.py -c POSCAR -r 0.2 -s 0.01 -o KPOINTS # -r: 相邻k点的距离(unit: ang), -s: 对称性检测的精度, -o: 指定输出名称
```
- 注意，该脚本来源于
> https://github.com/Chengcheng-Xiao/Tools/tree/master/VASP  
使用该脚本请关注 "https://github.com/Chengcheng-Xiao/"  
---
  
## 结构变换脚本  
1. MoveAlongAxis.py  
- 描述：沿着x,y,z轴整体移动原子  
- 环境：python, 需要安装numpy,ase,spglib  
``` shell
pip3 install numpy  
pip3 install ase 
pip3 install spglib
```
- 用法：  
``` shell
python MoveAlongAxis.py -c POSCAR -m 0 0 1 # 将POSCAR沿着z轴移动1A
python MoveAlongAxis.py -c POSCAR1 POSCAR2 -m 0 0 1 # 同时传入多个POSCAR
```
---
  
2. MirrorConversion.py  
- 描述：沿着指定的面镜像翻转所有原子
- 环境：python, 需要安装numpy,ase,spglib  
``` shell
pip3 install numpy  
pip3 install ase 
pip3 install spglib
```
- 用法：  
``` shell
python MirrorConversion.py -c POSCAR -x -y # 沿着(1 0 0)和(0 1 0)面镜像翻转
python MirrorConversion.py -c POSCAR1 POSCAR2 -x -y # 同时传入多个POSCAR
```
---
  
3. SortAtoms.py  
- 描述：沿着指定的轴对原子进行排序
- 环境：python, 需要安装numpy,ase,spglib,ordered_set  
``` shell
pip3 install numpy  
pip3 install ase 
pip3 install spglib
pip3 install ordered_set  
```
- 用法：  
``` shell
python SortAtoms.py -c POSCAR -x # 沿着x轴排序
python SortAtoms.py -c POSCAR1 POSCAR2 -x # 同时传入多个POSCAR
```
---
  
4. SuperCell.py  
- 描述：扩胞  
- 环境：python, 需要安装numpy,ase  
``` shell
pip3 install numpy  
pip3 install ase 
```
- 用法：  
``` shell
python SuperCell.py -c POSCAR -d --dim '2 2 2' # 扩222的胞 
```
- 注意，该脚本来源于
> https://github.com/Chengcheng-Xiao/Tools/tree/master/VASP  
使用该脚本请关注 "https://github.com/Chengcheng-Xiao/"  
---

5. ImposeSym.py  
- 描述：识别晶体对称性，并输出完善后的结构  
- 环境：python, 需要安装numpy,ase,spglib  
``` shell
pip3 install numpy  
pip3 install ase 
pip3 install spglib 
```
- 用法：  
``` shell
python ImposeSym.py -c POSCAR -s 0.1  # -s is the recongnize accuracy
```
- 注意，该脚本来源于
> https://github.com/Chengcheng-Xiao/Tools/tree/master/VASP  
使用该脚本请关注 "https://github.com/Chengcheng-Xiao/"  
---

6. ToFraction.py  
- 描述：将POSCAR转换成分数坐标  
- 环境：python, 需要安装numpy,ase  
``` shell
pip3 install numpy  
pip3 install ase 
```
- 用法：  
``` shell
python ToFraction.py -c POSCAR
```
---

7. ToCartesian.py  
- 描述：将POSCAR转换成笛卡尔坐标  
- 环境：python, 需要安装numpy,ase  
``` shell
pip3 install numpy  
pip3 install ase 
```
- 用法：  
``` shell
python ToCartesian.py -c POSCAR
```

# 后面会逐渐增加其他脚本的使用方法

