# VASP脚本手册
自用的VASP前处理和后处理脚本
转载时请注明出处
---
# 使用教程
1. Fixatoms.sh
- 描述：按照沿着z轴的顺序从底部开始固定胞原子
- 环境: Linux
- 用法：
``` shell
chmod u+x Fixatoms.sh
./Fixatoms.sh 16  # 固定POSCAR底部的16个原子
```

2. Fixatoms.py
- 描述：按照沿着z轴的顺序从底部开始固定胞原子
- 环境: python, 需要安装numpy
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

> 后面会逐渐增加其他脚本的使用方法

