import numpy as np
import ase
from ase.io import read
from ase.io import write
from ase.spacegroup import crystal
from ase.phonons import Phonons
import matplotlib.pyplot as plt
import matplotlib
import os
"""
描述：用ase调用vasp计算phonon

运行：
1. 更改参数设置
2. 在代码中修改vasp参数
3. 运行此文件

字体设置为Arial
如果系统没有Arial，下载Arial.ttf字体，使用以下命令获得matplotlib字体配置目录
>>>
import matplotlib.pyplot as plt
print(matplotlib.matplotlib_fname())  ---> /root/Software/anaconda3/envs/deepmd/lib/python3.10/site-packages/matplotlib/mpl-data/matplotlibrc
<<<
将Arial.ttf拷贝进./fonts/ttf/中
删除缓存 rm -rf ~/.cache/matplotlib
"""
# --------------参数设置-------------
plt.rcParams["font.sans-serif"] = ['Arial']
stru_path = "./cgN_unrelax.vasp" # 结构路径
calc_type = "VASP" # 计算器类型可以选"EMT","VASP","DP"
super_cell = (3,3,3) # 计算声子过程的扩胞大小
at_name = "attachment" #计算产生的文件存放的文件夹
phonon_color = "blue" # phonon颜色，可以是#111111
image_ratio = (7,5) # 图像比例
lwith = 2 # 线条粗细
bwith = 2 # 边框粗细
bias = -0.0005 # y轴上展示区域的偏移值，如果偏移值过大，图像在y轴上会过于拥挤
dos_fill_color = "grey" # dos填充颜色
dos_edge_color = "black" # dos边缘颜色
# --------------参数设置-------------

# 创建文件夹
if not os.path.exists(at_name) and not os.path.isdir(at_name):
    os.mkdir(at_name)

# 读取结构信息
atoms = read(stru_path)

# ------------VASP参数设置---------------
if calc_type.lower() == "vasp": 
    from ase.calculators.vasp import Vasp
    calc = Vasp(\
        directory="./workdir", # 计算路径
        txt="./out.log", # vasp输出文件名
        xc="pbe",
        npar=8,
        #lreal="auto",
        #isym=0,
        istart=0,
        icharg=2,
        prec="Accurate",
        ibrion=-1,
        #isif=3,
        #nsw=500,
        encut=520,
        nelm=300,
        nelmin=8,
        ediff=1E-6,
        #ediffg=ediffg,
        ismear=0,
        sigma=0.05,
        lwave=".F.",
        lcharg=".F.",
        kpts = [10,10,10], 
        gamma = "False",
        algo = "FAST",
        #ialgo=38,
        setups='recommended',\
        )
# ------------VASP参数设置---------------
elif calc_type.lower() == "dp": 
    from deepmd.calculator import DP
    calc = DP(model="graph.pb")
elif calc_type.lower() == "emt": 
    from ase.calculators.emt import EMT
    calc=EMT()
else:
    print("calc_type输入有误!!!")

# 用DP势绑定结构
# atoms.calc = calc

# 如果计算过声子，则不需要重新计算
os.system("if [ -d phonon_ ]; then cp -r ./phonon_ ./phonon; fi")

# 声子计算
ph = Phonons(atoms, calc, supercell=super_cell, delta=0.01)
ph.run()

# 保留声子计算过程中的数据
os.system("if [ -d phonon_ ]; then /bin/rm -r phonon_; fi; cp -r phonon phonon_")

# 读取力并组装动力学矩阵
ph.read(method='standard', acoustic=True) # frederiksen or standard
ph.clean()

# 获得K_path
path = atoms.cell.get_bravais_lattice().bandpath(npoints=100)  # path.path >>> "GXMGRX,MR"
# 获得声子带结构 
bs = ph.get_band_structure(path)

# 计算dos
dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)

# 绘图
def phonopy_plot(bs, ax, color="b", bwith=1, lwith=1, bais=bias, show=False, savepath=None): # savepath: 图像保存路径

    # 绘图数据
    energies = bs.energies  # energies.shape: 1,100,24 # 24条带
    xcoords, label_xcoords, orig_labels = map(list,bs.get_labels()) # xcoords.shape: 100
    # label_xcoords: [0.         0.83442036 1.66884072 2.84888931 4.29414777 5.47419635 5.47419635 6.30861671]
    # orig_labels: ['G', 'X', 'M', 'G', 'R', 'X', 'M', 'R']

    # 美化K路径符号
    def pretty(kpt):
        if kpt == 'G':
            kpt = r'$\Gamma$'
        elif len(kpt) == 2:
            kpt = kpt[0] + '$_' + kpt[1] + '$'
        return kpt

    labels = [pretty(name) for name in orig_labels] # ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X', 'M', 'R']

    # 这一步处理的目的是将重复的k_label写在一起，比如K和M的x轴重合了，写成K,M
    i = 1
    while i < len(labels):
        if label_xcoords[i - 1] == label_xcoords[i]:
            labels[i - 1] = labels[i - 1] + ',' + labels[i]
            labels.pop(i)
            label_xcoords.pop(i)
        else:
            i += 1

    # 绘制k点的竖线
    for x in label_xcoords[1:-1]:
        ax.axvline(x, color='0.5')

    # y label
    ylabel = 'Energies [eV]'

    # x tick
    ax.set_xticks(label_xcoords)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)

    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    # 在能量参考处画一条虚线
    ax.axhline(bs.reference, color='k', ls=':')

    # 图像范围
    emin= dos.get_energies().min()-bais
    emax = dos.get_energies().max()+bais
    ax.axis(xmin=0, xmax=xcoords[-1], ymin=emin, ymax=emax)

    for spin, e_kn in enumerate(energies): # e_kn.shape: 100,24
        ax.plot(xcoords, e_kn[:, 0], color=color, lw=lwith)
        # 绘制每一条能带
        for e_k in e_kn.T:
            ax.plot(xcoords, e_k, color=color, lw=lwith)
    if show==True:
        plt.show()
    if savepath:
        plt.savefig(savepath)

image_ratio = ((image_ratio[0]/image_ratio[1]) * 5, 5)
fig = plt.figure(1, figsize=image_ratio, dpi=100)
ax = fig.add_axes([.12, .07, .63, .85]) # 图像范围
phonopy_plot(bs, ax, color=phonon_color, bwith=bwith, lwith=lwith, show=False, savepath=None)

# 绘制dos
dosax = fig.add_axes([.75, .07, .17, .85])
dosax.spines['bottom'].set_linewidth(bwith)
dosax.spines['left'].set_linewidth(bwith)
dosax.spines['top'].set_linewidth(bwith)
dosax.spines['right'].set_linewidth(bwith)

dosax.fill_between(x=dos.get_weights(), y1=dos.get_energies(), y2=0, color=dos_fill_color,edgecolor=dos_edge_color, lw=lwith)

dosax.set_ylim(dos.get_energies().min()-bias, dos.get_energies().max()+bias)
dosax.set_xlabel("DOS", fontsize=10)
dosax.set_yticks([])
dosax.set_xticks([])

# 在能量参考处处画一条虚线
dosax.axhline(bs.reference, color='k', ls=':',lw=2)

fig.savefig(f"{at_name}/8_phonon_.png")

