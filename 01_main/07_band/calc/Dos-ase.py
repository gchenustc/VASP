from ase.io import read
from ase.calculators.vasp import Vasp
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
import os
"""
描述：用ase调用vasp计算dos

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
plt.rcParams["font.sans-serif"] = ['Arial']

#  --------------- 参数设置 --------------- 
stru_path = "./POSCAR" # 结构路径
relax = True # 是否进行结构优化
vasp_calc  = False# 是否调用vasp计算，如果是False则在vasp_calc_dir目录中获取vasp已经计算好的结果进行绘图
vasp_calc_dir = './12_dos' # vasp工作目录路径
image_ratio = (6,5) # 图像比例
color = "blue" # 能带颜色，可以是#111111
lwith = 2 # 线条粗细
bwith = 2 # 边框粗细
bias = 0.2
at_name = "attachment" # 计算产生的文件存放的文件夹
#  --------------- 参数设置 --------------- 

# 创建文件夹
if not os.path.exists(at_name) and not os.path.isdir(at_name):
    os.mkdir(at_name)

if not vasp_calc:
    calc = Vasp(restart=True, directory=vasp_calc_dir)

else:
    print("iter0>>>>>>")
    # 设置vasp结构优化
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
    atoms = read("./CONTCAR") if relax else read(stru_path)
    # 设置vasp自洽参数
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
        lorbit=11,
        setups='recommended',\
        )
    # 将原子和计算器绑定
    atoms.calc = calc
    # 计算
    energy = atoms.get_potential_energy()
    print("energy:",energy)

    print("iter2>>>>>>")
    atoms = read("./CONTCAR")
    # 设置vasp非自洽参数
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
        kpts = [10, 10, 10],
        gamma = False,
        algo = "FAST",
        #ialgo=38,
        #ispin=2,
        lorbit=10,
        setups='recommended',\
        )
    # 将原子和计算器绑定
    atoms.calc = calc
    # 计算
    energy = atoms.get_potential_energy()
    print("energy:",energy)

# DOS
dos = DOS(calc, width=0.2)
d = dos.get_dos()
e = dos.get_energies()

# plot
image_ratio = ((image_ratio[0]/image_ratio[1]) * 5, 5)
fig = plt.figure(1, figsize=image_ratio, dpi=100)
ax = fig.add_subplot(111)
ax.axis(xmin=e.min(), xmax=e.max(), ymin=d.min()-0.01, ymax=d.max()+bias)
ax.plot(e, d, color=color, lw=lwith)
ax.set_xlabel('energy [eV]')
ax.set_ylabel('DOS')
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
plt.savefig(f"{at_name}/10_dos.png")
