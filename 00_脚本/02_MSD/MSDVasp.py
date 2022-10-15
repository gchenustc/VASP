"""
0.需要POSCAR和XDATCAR文件（晶格常数从XDATCAR获取，由于 ISIF=2，晶胞不改变。用 XDATCAR 种的第一步作为初始结构）
1.不需要 INCAR，因为步长是手动输入的
2.需要安装python，默认的python命令是python，如果需要更改，请至msd_vasp.sh文件将"python python_msd.py"这一行更改
3.需要numpy包 --> pip install numpy
4.如果需要导出图片，需要matplotlib,pandas包 --> pip install pandas, pip install matplotlib
并且需要在Linux中创建文件 ~/.config/matplotlib/matplotlibrc，（其中，~/.config/matplotlib/是配置文件matplotlibrc的路径）。在下面一行添加
backend : Agg

注意：
只适用于计算立方格子的MSD
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib as mpl

#--------------------define function--------------------

def fra2car(pos,cons):
    """将分数坐标转换为笛卡尔坐标，传入的参数分别为分数坐标和晶格常数"""
    return np.dot(pos,cons.T)

def calc_msd(file1,file2,demision=3,pmsd_list=None):
    """计算msd，
    demision指的是维度，1指的是算每一个维度上的msd,返回的是元组，demision只能取1或者3
    pmsd_list指的是传入包含每种原子个数的列表，比如[3,9]，用于计算分msd"""
    
    assert file1.shape == file2.shape
    
    if pmsd_list:
        result=[]
        for i,j in zip(np.vsplit(file1,np.cumsum(pmsd_list)[:-1]),np.vsplit(file2,np.cumsum(pmsd_list)[:-1])):
            result.append(calc_msd(i,j))
        return result
    
    if demision == 3:
        return np.square(file1-file2).sum(axis=1).mean()
    elif demision == 1:
        return np.square(file1-file2).mean(axis=0)

def judge_period_atoms(before_pos,after_pos,cons,record):
    """原子位置的周期性判断,cons指的是晶格常数"""
    for i in range(3):
        pos_diff = after_pos[:,i] - before_pos[:,i]
        judge_condition =  0.5 * cons[i,i]
        record[:,i] += np.where(pos_diff < -judge_condition,1,0)
        record[:,i] += np.where(pos_diff > judge_condition,-1,0)

def calibration_pos(pos,cons,record):
    """根据周期性边界条件的表格校准位置"""
    for i in range(3): #3d
        pos[:,i] += record[:,i]*cons[i,i]
    return pos

def plot():
    #配置
    mpl.use('Agg')
    plt.rcParams['font.sans-serif']=['Arial']
    formatter = ticker.FormatStrFormatter("%.2e")
    msd = pd.read_csv("msd.csv",sep=" ",index_col="time(fs)")
    #绘制图1
    fig = plt.figure(dpi=400,figsize=(10,6.8))
    ax = fig.add_subplot(111)
    ax.plot(msd.index,msd.iloc[:,0:4],aa=True)
    ax.legend(msd.columns,frameon=False,loc="upper center",ncol=4)
    ax.set_title("MSD",fontsize=15)
    ax.set_xlabel(msd.index.name)
    ax.set_ylabel("MSD(A^2)")
    ax.set_xlim(0,msd.index[-1])
    ax.yaxis.set_major_formatter(formatter)
    plt.savefig("msd.png")
    plt.show()
    plt.close()
    #绘制图2
    fig1,ax1 = plt.subplots(1,1,dpi=400,figsize=(10,6.8))
    ax1.plot(msd.index,msd.iloc[:,4:],aa=True)
    ax1.legend(msd.iloc[:,4:].columns,frameon=False,loc="upper center",ncol=len(msd.iloc[:,4:].columns))
    ax1.set_title("MSD",fontsize=15)
    ax1.set_xlabel(msd.index.name)
    ax1.set_ylabel("MSD(A^2)")
    ax1.set_xlim(0,msd.index[-1])
    ax1.yaxis.set_major_formatter(formatter)
    plt.savefig("msd-element.png")
    plt.show()
    plt.close()
    
#--------------------end define function--------------------

def main():
    # 导入文件
    constant = np.loadtxt("cons.csv",dtype='f8') #晶格常数 - 来自于 XDATCAR 的第一步的晶格常数，由于 MD 计算是固定胞的，所有结构胞都一样
    
    atoms_info = np.loadtxt("atoms_info.csv",dtype=str)  #原子种类信息
    #atoms_info = np.array(['N','192'])  # 模拟胞内只有一种原子的情况
    
    # 如果胞内只有一种原子，atoms_info就会如 ['N','192'] 这样的一维形式，把它换成二维 [['N'],['192']]
    if atoms_info.ndim == 1:
        atoms_info = atoms_info[:,None]
        
    atoms_num = sum([int(i) for i in atoms_info[1]])    # 一个晶胞内的原子数
    
    init_pos_fra = np.loadtxt("init_pos.csv") # 分数坐标的POSCAR位置坐标
    xdatcar_fra =  np.loadtxt("xdatcar.csv") # 分数坐标的xdatcar位置坐标
    steps,step_len =  np.loadtxt("step.csv")  #步数和步长，步长是用户输入的，步数是从XDATCAR 导入的
    
    # 将 分数坐标 转换为 笛卡尔坐标
    init_pos = fra2car(init_pos_fra,constant)
    xdatcar = np.vsplit(fra2car(xdatcar_fra,constant),steps)
    init_pos = xdatcar[0] # 用xdatcar的第1步作为初始位置
    
    #print(constant)
    #print(atoms_info)
    #print(atoms_num)
    #print(init_pos_fra)
    #print(xdatcar_fra)
    
    # 设置原子信息字典
    atoms_info_dict={}  #字典中存放原子信息 eg {N:192,H:8}
    for i in range(atoms_info.shape[1]):
        atoms_info_dict.setdefault(atoms_info[0,i],int(atoms_info[1,i]))
    
    #print(atoms_info_dict)
    
    
    # 定义一些变量
    msd_list = ['time(fs)','MSD-total','MSD-x','MSD-y','MSD-z']  #创建msd输出文件的头行
    for i in atoms_info_dict.keys():
        msd_list.append('MSD-'+i)
    #print(msd_list)
    
    over_boundary_record = np.zeros((atoms_num,3))  #记录原子超出边界的次数，纵坐标对应每个原子，横坐标使对应每个原子的 x,y,z 坐标，某个原子的某个轴超出边界，对应位置+1或者-1
    #print(over_boundary_record)
    
    
    # ------------ 信息导入和处理完毕，下面是主要部分 --------------
    for step in range(len(xdatcar)):
        pos = xdatcar[step]
        if step != 0:
            judge_period_atoms(xdatcar[step-1],pos,constant,over_boundary_record)
        pos=calibration_pos(pos,constant,over_boundary_record)
        over_boundary_record = np.zeros((atoms_num,3)) #初始化
        #print(step+1)
        #print(pos)
        #print(over_boundary_record)
        #print()
        
        msd_list_append = [] #输出文件的一行
        msd_list_append.append((step+1)*step_len) #输出当前时间
        msd_list_append.append(calc_msd(pos,init_pos,3)) #输出总msd
        msd_list_append.extend(calc_msd(pos,init_pos,1)) #输出pmsd
        msd_list_append.extend(calc_msd(pos,init_pos,3,list(atoms_info_dict.values()))) #输出pmsd
        msd_list  = np.vstack([msd_list,msd_list_append])
        #print(step+1)
        #print(pos)
        #print(msd_list)

    np.savetxt("msd.csv",msd_list,fmt="%s")
    plot()  #绘图

if __name__ == '__main__':
    main()
