import numpy as np
import ase
from ase.io import read
from ase.io import write
from ase.spacegroup import crystal
from ase.calculators.vasp import Vasp
import matplotlib
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
"""
描述：已经用vasp计算好band和dos，利用此脚本绘图

运行：
1. 更改参数设置
2. 运行此文件

字体设置为Arial
如果系统没有Arial，下载Arial.ttf字体，使用以下命令获得matplotlib字体配置目录
>>>
import matplotlib.pyplot as plt
print(matplotlib.matplotlib_fname())  ---> /root/Software/anaconda3/envs/deepmd/lib/python3.10/site-packages/matplotlib/mpl-data/matplotlibrc
<<<
将Arial.ttf拷贝进./fonts/ttf/中
删除缓存 rm -rf ~/.cache/matplotlib
"""
#  --------------- 参数设置 --------------- 
plt.rcParams["font.sans-serif"] = ['Arial']
band_calc_dir = "11_band"  # 计算能带的文件夹
dos_calc_dir = "12_dos"  #　计算DOS的文件夹
emin=-6 # y轴最小值
emax=15 # y轴最大值
image_ratio = (7,5) # 图像比例
color = "blue" # 能带颜色，可以是#111111
lwith = 2 # 线条粗细
bwith = 2 # 边框粗细
dos_fill_color = "grey"
dos_edge_color = "black"
at_name = "attachment" # 计算产生的文件存放的文件夹
bias=0.5
#  --------------- 参数设置 --------------- 

# 创建文件夹
if not os.path.exists(at_name) and not os.path.isdir(at_name):
    os.mkdir(at_name)

calc1 = Vasp(restart=True, directory=band_calc_dir)
calc2 = Vasp(restart=True, directory=dos_calc_dir)

# 获得能带
bs = calc1.band_structure()
# 绘图 - 自带的函数绘制，不好自定义
# bs.plot(emin=0, emax=20, filename='band.png')

# 绘图
# 图像位置
def band_plot(bs, ax, emin=None, emax=None, color="b", bwith=1, lwith=1, show=False, savepath=None): # savepath: 图像保存路径

    # 绘图数据
    energies = bs.energies - bs.reference  # energies.shape: 1,100,24 # 24条带
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
    ax.axhline(0, color='k', ls=':')

    
    if emin or emax:
        ax.axis(xmin=0, xmax=xcoords[-1], ymin=emin, ymax=emax)
    else:
        ax.axis(xmin=0, xmax=xcoords[-1], ymin=energies.min()-bias, ymax=energies.max()+bias)

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
ax1 = fig.add_axes([.12, .07, .63, .85]) # 图像范围

band_plot(bs, ax1, emin=emin, emax=emax, color=color,  bwith=bwith, lwith=lwith, show=False)
 
# dos
dos = DOS(calc2, width=0.2)
d = dos.get_dos()
e = dos.get_energies()

ax2 = fig.add_axes([.75, .07, .17, .85])
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)

ax2.fill_between(x=d, y1=e, y2=0, color=dos_fill_color,edgecolor=dos_edge_color, lw=lwith)

ax2.set_ylim(emin, emax)
ax2.set_xlabel("DOS", fontsize=10)
ax2.set_yticks([])
ax2.set_xticks([])

# 在能量参考处处画一条虚线
ax2.axhline(0, color='k', ls=':',lw=2)
fig.savefig(f"{at_name}/11_band_dos.png")