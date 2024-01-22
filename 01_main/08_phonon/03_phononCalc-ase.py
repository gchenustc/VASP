"""
Description: Calculating phonons using VASP with ASE.

Usage:
1. Modify the parameter settings.
2. Modify the VASP parameters in the code.
3. Run this file.

Note:
1. This program needs to be run on a Linux system.
2. ASE needs to be installed, and the VASP environment required by ASE (including the location of the potential function and the VASP command) or the environment for DP potential needs to be set up.
3. The font is set to Arial. If Arial is not available on your system, download the Arial.ttf font and use the following command to obtain the matplotlib font configuration directory:
   >>>
   import matplotlib.pyplot as plt
   print(matplotlib.matplotlib_fname())  ---> /root/Software/anaconda3/envs/deepmd/lib/python3.10/site-packages/matplotlib/mpl-data/matplotlibrc
   <<<
   Copy Arial.ttf to ./fonts/ttf/.
   Delete the cache: rm -rf ~/.cache/matplotlib
"""

from ase.io import read
from ase.phonons import Phonons
import matplotlib.pyplot as plt
import os

# --------------Parameter Settings-------------
plt.rcParams["font.sans-serif"] = ['Arial']
stru_path = "./POSCAR"  # Structure path
calc_type = "VASP"  # Calculator type, can be "EMT", "VASP", or "DP"
super_cell = (3, 2, 1)  # Supercell size for phonon calculation
at_name = "attachment"  # Folder to store the generated files

# --------------Graph Parameter Settings-------------
phonon_color = "blue"  # Phonon color, can be specified as "#111111"
image_ratio = (7, 5)  # Image ratio
lwith = 2  # Line width
bwith = 2  # Border width
bias = -0.0005  # Offset value for the y-axis display, if the offset value is too large, the image will be crowded on the y-axis
y_label_size = 12  # Size of the y-axis label
dos_fill_color = "grey"  # DOS fill color
dos_edge_color = "black"  # DOS edge color

# --------------Parameter Settings-------------

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
        setups={"N":"_h", "F":"_h"},
        npar=8, 
        lreal=False,
        #isym=0,
        istart=0,
        icharg=2,
        prec="Accurate",
        #addgrid=True,  # 增加精度
        ibrion=-1,
        #nfree=2,
        #isif=3,
        #nsw=500,
        encut=1000,
        nelm=300,
        nelmin=8,
        ediff=1E-8,
        pstress=2000,   # change!!!!!!
        #ediffg=-1e-3,
        ismear=0,   # -5
        sigma=0.05,
        lwave=False,
        lcharg=False,
        kspacing=0.15, 
        kgamma=False,
        #algo="FAST",
        ialgo=38,
        ivdw=11,
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
path = atoms.cell.get_bravais_lattice().bandpath(npoints=200)  # path.path >>> "GXMGRX,MR"
# 获得声子带结构 
bs = ph.get_band_structure(path)

# 计算dos
dos = ph.get_dos(kpts=(40, 40, 40)).sample_grid(npts=200, width=1e-3)

# 绘图
def phonopy_plot(bs, ax, color="b", bwith=1, lwith=1, bais=bias, show=False, savepath=None): # savepath: 图像保存路径

    # 绘图数据
    frenquency = bs.energies/0.0041  # energies.shape: 1,100,24 # 24条带
    xcoords, label_xcoords, orig_labels = map(list,bs.get_labels()) # xcoords.shape: 100
    # label_xcoords: [0.         0.83442036 1.66884072 2.84888931 4.29414777 5.47419635 5.47419635 6.30861671]
    # orig_labels: ['G', 'X', 'M', 'G', 'R', 'X', 'M', 'R']

    # 美化K路径符号
    def pretty(kpt):
        if kpt == 'G':
            kpt = r'$\Gamma$'
        elif len(kpt) == 2:
            kpt = f'{kpt[0]}$_{kpt[1]}$'
        return kpt

    labels = [pretty(name) for name in orig_labels] # ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X', 'M', 'R']

    # 这一步处理的目的是将重复的k_label写在一起，比如K和M的x轴重合了，写成K,M
    i = 1
    while i < len(labels):
        if label_xcoords[i - 1] == label_xcoords[i]:
            labels[i - 1] = f'{labels[i - 1]},{labels[i]}'
            labels.pop(i)
            label_xcoords.pop(i)
        else:
            i += 1

    # 绘制k点的竖线
    for x in label_xcoords[1:-1]:
        ax.axvline(x, color='0.5')

    # y label
    ylabel = 'Frenquency (THZ)'

    # x tick
    ax.set_xticks(label_xcoords)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel, fontsize=y_label_size)

    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    # 在能量参考处画一条虚线
    ax.axhline(bs.reference, color='k', ls=':')

    # 图像范围
    emin= (dos.get_energies().min())/0.0041-bais
    emax = (dos.get_energies().max())/0.0041+bais
    ax.axis(xmin=0, xmax=xcoords[-1], ymin=emin, ymax=emax)

    for spin, e_kn in enumerate(frenquency): # e_kn.shape: 100,24
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
dosax.set_xlabel("DOS", fontsize=12)
dosax.set_yticks([])
dosax.set_xticks([])

# 在能量参考处处画一条虚线
dosax.axhline(bs.reference, color='k', ls=':',lw=2)

# 保存文件
fig.savefig(f"{at_name}/phonon.png")