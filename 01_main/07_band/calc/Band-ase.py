from ase.io import read
from ase.calculators.vasp import Vasp
import matplotlib.pyplot as plt
import os
"""
描述：用ase调用vasp计算band

运行：
1. 更改参数设置
2. 在代码中修改vasp参数
3. 运行此文件

字体设置：
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
plt.rcParams["font.sans-serif"] = ['Arial'] # 字体设置
stru_path = "./POSCAR" # 结构路径
relax = False # 是否进行结构优化
vasp_calc = False # 是否调用vasp计算，如果是False则在vasp_calc_dir目录中获取vasp已经计算好的结果进行绘图
vasp_calc_dir = './11_band' # vasp工作目录路径
emin=-6
emax=15
image_ratio = (6,5) # 图像比例
color = "blue" # 能带颜色，可以是#111111
lwith = 2 # 线条粗细
bwith = 2 # 边框粗细
at_name = "attachment"
bias=0.5
#  --------------- 参数设置 ---------------

# 创建文件夹
if not os.path.exists(at_name) and not os.path.isdir(at_name):
    os.mkdir(at_name)

if not vasp_calc:
      # 读取vasp自洽的结果
      calc = Vasp(restart=True, directory=vasp_calc_dir)
else:
    print("iter0>>>>>>")

# --------------- 设置vasp结构优化 ---------------
    if relax:
        atoms = read(stru_path)
        calc = Vasp(\
            directory="./", # 计算路径
            txt="./out.log", # vasp输出文件名
            xc="pbe",
            npar=8,
            #lreal="auto",
            #isym=0,
            istart=0,
            icharg=2,
            prec="Accurate",
            ibrion=2,
            isif=2,
            nsw=500,
            encut=520,
            nelm=100,
            nelmin=8,
            ediff=1E-6,
            ediffg=1E-3,
            ismear=0,
            sigma=0.05,
            lwave=".F.",
            lcharg=".F.",
            kpts = [10, 10, 10],
            gamma = False,
            algo = "FAST",
            #ialgo=38,
            #ispin=2,
            #lorbit=11,
            setups='recommended',\
            )
        # 将原子和计算器绑定
        atoms.calc = calc
        # 计算
        energy = atoms.get_potential_energy()
        print("energy:",energy)

    print("iter1>>>>>>")
    if relax:
        atoms = read("./CONTCAR")
    else:
        atoms = read(stru_path) 

# --------------- 设置vasp自洽计算 ---------------
    calc = Vasp(\
        directory="./", # 计算路径
        txt="./out.log", # vasp输出文件名
        xc="pbe",
        npar=8,
        #lreal="auto",
        #isym=0,
        istart=0,
        icharg=2,
        prec="Accurate",
        ibrion=-1,
        #isif=2,
        #nsw=500,
        encut=520,
        nelm=100,
        nelmin=8,
        ediff=1E-6,
        ediffg=1E-3,
        ismear=0,
        sigma=0.05,
        lwave=".T.",
        lcharg=".T.",
        kpts = [10, 10, 10],
        gamma = False,
        algo = "FAST",
        #ialgo=38,
        #ispin=2,
        #lorbit=11,
        setups='recommended',\
        )
    # 将原子和计算器绑定
    atoms.calc = calc
    # 计算
    energy = atoms.get_potential_energy()
    print("energy:",energy)

    nbands = int(os.popen("echo `awk '/NBANDS/{print $NF}' OUTCAR`").read().strip())
    print("nbands",nbands)

    print("iter2>>>>>>")
    atoms = read("./CONTCAR")
    # 获得高对称点路径
    bs = atoms.cell.get_bravais_lattice().bandpath(npoints=300)

# --------------- 设置vasp非自洽计算 ---------------
    calc = Vasp(\
        directory="./", # 计算路径
        txt="./out.log", # vasp输出文件名
        xc="pbe",
        npar=8,
        #lreal="auto",
        #isym=0,
        istart=1,
        icharg=11,
        prec="Accurate",
        ibrion=-1,
        #isif=2,
        #nsw=500,
        encut=520,
        nelm=100,
        nelmin=8,
        ediff=1E-6,
        ediffg=1E-3,
        ismear=0,
        sigma=0.05,
        lwave=".T.",
        lcharg=".T.",
        kpts = {"path":bs.path,"npoints":100},
        gamma = False,
        algo = "FAST",
        #ialgo=38,
        #ispin=2,
        nbands=nbands,
        setups='recommended',\
        )
    # 将原子和计算器绑定
    atoms.calc = calc
    # 计算
    energy = atoms.get_potential_energy()
    print("energy:",energy)

# 获得能带
bs = calc.band_structure()

# 绘图 - 自带的函数绘制，不好自定义
# bs.plot(emin=0, emax=20, filename='band.png')

# 绘图
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
ax = fig.add_subplot(111)

band_plot(bs, ax, emin=emin, emax=emax, color=color,  bwith=bwith, lwith=lwith, show=False, savepath=f"./{at_name}/9_band_.png")
