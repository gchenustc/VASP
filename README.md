# VASP脚本手册
---
如果有bug，欢迎加我wechat反馈: c1337375425

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

3. MSD.py
- 描述：计算分子动力学的MSD, 需要准备XDATCAR, 输出msd.png, msd-element.png, msd.out, msd.out文本的形式可以在origin绘图
- 环境：python
python需安装numpy，matplotlib，pandas和ase
``` sheLL
pip install numpy
pip install matplotlib
pip install pands
pip install ase
```
- 用法：
```
python MSD.py -p 0.5 # -p 指定步长
```
- 其他:
- MSD_tetra.py 与之用法一样，但只支持四方晶体
- MSD_Liujiapeng.py 该脚本来源于:
> https://github.com/jiapeng-liu/Mean-Squared-Displacment_msd_for_vasp
---

4. chgdiff.py
- 描述：计算差分电荷
- 环境：python, 需要安装numpy,ase
``` shell
pip3 install numpy
pip3 install ase
```
- 用法：
``` shell
python chgdiff.py CHGCAR1 CHGCAR2  #输出文件名为 CHGDIFF.
```
- 注意，该脚本来源于
> https://github.com/Chengcheng-Xiao/Tools/tree/master/VASP
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
---

6. CalcCohp.py
- 描述：计算cohp
- 环境：python, 需要安装numpy,ase
``` shell
pip3 install numpy
pip3 install ase
```
- 用法：
``` shell
python current_file.py -l lobster  # -l 后指定lobster的执行命令
python current_file.py -l lobster -p -a 2 4 5 7  # 出现 -p 则表示计算指定原子的cohp(否则计算所有原子的cohp), -a 后为指定原子的序号，比如示例将计算序号为2和4以及5和7原子的cohp
```
- 注意，用ase调用vasp计算需要配置好vasp的执行环境和POTCAR的目录位置，可以看官网教程
> https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#introduction
---

7. PlotCohp.py
- 描述：用CalcCohp.py计算cohp后，用此脚本批量绘图
- 环境：python, 需要安装numpy,matplotlib
``` shell
pip3 install numpy
pip3 install matplotlib
```
- 用法：
``` shell
python PlotCohp.py
```
---

8. AveBondLen.py
- 描述：按照成键种类计算其平均键长，当然也给出所有键的平均键长
- 环境：python, 需要安装numpy,ase
``` shell
pip3 install numpy
pip3 install ase
```
- 用法：
``` shell
python AveBondLen.py
```
---

9. Band-Ase.py
- 描述：用ASE调用VASP去计算能带，并绘制图形  
- 环境：python, 需要安装numpy,ase,matplotlib   
``` shell
pip3 install numpy
pip3 install ase
pip3 install matplotlib
```
- 用法：
1. 在代码中调整参数 
2. 在代码中修改vasp的INCAR设置  
3. 运行此文件  
``` shell
python Band-Ase.py
```
---

10. Dos-Ase.py
- 描述：用ase调用vasp计算dos，并绘制图形  
- 环境：python, 需要安装numpy,ase,matplotlib  
``` shell
pip3 install numpy
pip3 install ase
pip3 install matplotlib
```
- 用法： 
1. 在代码中调整参数 
2. 在代码中修改vasp的INCAR设置  
3. 运行此文件 
``` shell
python Dos-Ase.py
```
---

11. Phonon-Ase.py
- 描述：用ase调用vasp计算phonon，并绘制图形   
- 环境：python, 需要安装numpy,ase,matplotlib  
``` shell
pip3 install numpy
pip3 install ase
pip3 install matplotlib
```
- 用法：
1. 在代码中调整参数 
2. 在代码中修改vasp的INCAR设置  
3. 运行此文件  

- 注意： 
1. 需要在linux系统下运行该程序  
2. 需要安装ase，并且有ase需要的vasp环境（包括势函数位置和vasp运行的命令）或者DP势的环境（关于如何用ase调用vasp请看官网教程） 
3. 默认将图片的字体设置为Arial  
如果系统没有Arial，请下载Arial.ttf字体，否则用默认字体，使用以下命令获得matplotlib字体配置目录  
```python
import matplotlib.pyplot as plt
print(matplotlib.matplotlib_fname())
```
输出  
``` shell
/root/Software/anaconda3/envs/deepmd/lib/python3.10/site-packages/matplotlib/mpl-data/matplotlibrc
```
将Arial.ttf拷贝进./fonts/ttf/中  
删除缓存 rm -rf ~/.cache/matplotlib后运行脚本  
``` shell
python Phonon-Ase.py
```
---

12. StruEquiTest.py  
- 描述：检测两个结构是否为等价，支持多种文件格式     
- 环境：python, 需要安装numpy,ase,pymatgen  
``` shell
pip3 install numpy
pip3 install ase
pip3 install pymatgen
```
- 用法：
``` shell
python *.py POSCAR1 POSCAR2 ...  # POSCAR*为结构文件的路径
```
---

13. doping.py  
- 描述：随机替换掺杂，需要提供模板结构，并更改 “tem_stru_path” 这个变量为结构文件名，其他参数请看文件的末尾。  
在“old_strus_dir_name”中的>子文件夹<中的结构为排除的结构，脚本生成的结构会排除与其等价的结构。可以不提供该文件夹，但此变量必须存在。  
不要将排除的结构直接放在old_strus_dir_name下 ，需要在此文件夹中创建子文件夹，然后再将排除的结构放在子文件夹中，比如old_strus_dir_name=old_strus  
```shell
tree old_strus
```
>>>
old_strus/a/a.vasp  
old_strus/a/b.vasp  
old_strus/b/a.vasp  
old_strus/b/b.vasp  
<<<

- 环境：python, 需要安装numpy,ase,pymatgen   
``` shell
pip3 install numpy
pip3 install ase
pip3 install pymatgen
```
用法： 
与general_fun放在同一个文件夹下，直接运行该文件  
---

14. CalcEnergyDensity.py 
- 描述：计算一个结构的能量密度，需要给定所有基态反应物单质的能量（平均到单个原子上），修改energy_reaction参数。比如如果要计算C3N4的含能量，C基态结构的能量为-3eV/atom, N的基态能量为-5eV/atom，则energy_reaction={"C":-3, "N":-5}，其他变量不需要修改
- 环境：python, 需要安装pymatgen   
``` shell
pip3 install pymatgen
```
用法：
和POCAR和OSZICAR在同一个文件夹下，然后直接运行该文件
---

15. PosEquiText.py  
- 描述：检测POSCAR结构中的两个原子位点否为等价     
- 环境：python, 需要安装numpy,ase,pymatgen  
``` shell
pip3 install numpy
pip3 install ase
pip3 install pymatgen
```
- 用法：
``` shell
python *.py  # 需要在脚本内部设定参数 1. vasp结构文件名，2. 指定检测的某种元素类别
```
---

16. Calypso_sta.py & Magus_sta.py  
- 描述：请看脚本注释
---

## 结构变换脚本
1. MoveAlongAxis.py
- 描述：沿着x,y,z轴（笛卡尔坐标）整体移动原子
- 环境：python, 需要安装numpy,ase,spglib
``` shell
pip3 install numpy
pip3 install ase
pip3 install spglib
```
- 用法：
``` shell
python *.py -c POSCAR1 POSCAR2 -m 0 0 1 # 将POSCAR1和POSCAR2沿着z轴（Cartesian坐标）移动1A，可以不传入-c，默认文件名为POSCAR
python *.py -c POSCAR1 POSCAR2 -d direct -m 0 0 1 -v # -d指移动mode，沿着c轴（分数坐标）移动
python *.py -c POSCAR1 POSCAR2 -d direct -p -m 0 0 0.1 -v # -p打开后，按照胞长的比例移动，-p参数需要-d direct
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
---





