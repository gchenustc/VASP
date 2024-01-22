"""
Created on 2022-11-4
author: gchen
Contact: gchenn@mail.ustc.edu.cn
introductions:
must install numpy, ase, matplotlib for python3
```
pip install numpy,ase,matplotlib
```
call script:
```
python file.py -p 0.5 # -p is lenth per step
```
"""

import numpy as np
import pandas as pd
import ase
import ase.io
import argparse
from ase.visualize import view
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib as mpl
import time

# 绘图函数
def plot():
    #配置
    mpl.use('TkAgg')
    plt.rcParams['font.sans-serif']=['Arial']
    formatter = ticker.FormatStrFormatter("%.2e")
    msd = pd.read_csv("msd.out",sep=" ",index_col="time(fs)")
    #绘制图1
    fig = plt.figure(dpi=60,figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.plot(msd.index,msd.iloc[:,0:4],aa=True)
    ax.legend(msd.columns,frameon=False,loc="upper center",ncol=4)
    ax.set_title("MSD",fontsize=15)
    ax.set_xlabel(msd.index.name)
    ax.set_ylabel("MSD(A^2)")
    ax.set_xlim(0,msd.index[-1])
    ax.yaxis.set_major_formatter(formatter)
    plt.savefig("msd.png")
    plt.close()
    #绘制图2
    fig1,ax1 = plt.subplots(1,1,dpi=60,figsize=(10,8))
    ax1.plot(msd.index,msd.iloc[:,4:],aa=True)
    ax1.legend(msd.iloc[:,4:].columns,frameon=False,loc="upper center",ncol=len(msd.iloc[:,4:].columns))
    ax1.set_title("MSD",fontsize=15)
    ax1.set_xlabel(msd.index.name)
    ax1.set_ylabel("MSD(A^2)")
    ax1.set_xlim(0,msd.index[-1])
    ax1.yaxis.set_major_formatter(formatter)
    plt.savefig("msd-element.png")
    plt.close()

# 判断周期性
def judge_period_atoms(before_pos,after_pos,cons,record):
    """原子位置的周期性判断,cons指的是晶格常数"""
    for i in range(3):
        pos_diff = after_pos[:,i] - before_pos[:,i]
        judge_condition =  0.3 * cons[i,i]
        #if np.any(pos_diff < -judge_condition):print("1111")
        #if np.any(pos_diff > judge_condition):print("1111")
        record[:,i] += np.where(pos_diff < -judge_condition,1,0)
        record[:,i] += np.where(pos_diff > judge_condition,-1,0)

def calibration_pos(pos,cons,record):
    """根据周期性边界条件的表格校准位置"""
    for i in range(3): #3d
        pos[:,i] = pos[:,i] + record[:,i]*cons[i,i]
    return pos

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

if __name__ == "__main__":
    start_time = time.perf_counter()

    parser = argparse.ArgumentParser(description='give step lenth')
    parser.add_argument("-p","--potim", action="store", type=float, required=True, help="give potim for incar")
    options = parser.parse_args()

    from ase.io.vasp import read_vasp_xdatcar
    xdatcar = read_vasp_xdatcar("XDATCAR",index=0)
    xdatcar_0 = xdatcar[0]

    # 原子种类 -- eg. ['H', 'N']
    species = list(set(xdatcar_0.get_chemical_symbols()))
    # 每个原子种类对应的原子数目 -- eg. [12, 192]
    n_atoms = [len([atom for atom in xdatcar_0 if atom.symbol==specie]) for specie in species]
    n_atoms_total = sum(n_atoms)
    # 步长
    optim = options.potim
    # 第一个结构的坐标 - Cartesian坐标
    init_pos = xdatcar_0.get_positions()
    # 字典中存放原子信息 eg {"N":192,"H":8}
    atoms_info_dict={}
    for index,specie in enumerate(species):
        atoms_info_dict.setdefault(specie, n_atoms[index])

    # 输出
    msd_list = ['time(fs)','MSD-total','MSD-x','MSD-y','MSD-z']  #创建msd输出文件的头行
    for i in atoms_info_dict.keys():
        msd_list.append('MSD-'+i)

    #记录原子超出边界的次数，纵坐标对应每个原子，横坐标使对应每个原子的 x,y,z 坐标，某个原子的某个轴超出边界，对应位置+1或者-1
    over_boundary_record = np.zeros((n_atoms_total,3)) 

    pos_bak=None
    for n_step,atoms in enumerate(xdatcar):
        pos = atoms.get_positions()
        
        if n_step != 0:
            # 前一步的pos
            judge_period_atoms(pos_bak, pos, xdatcar_0.cell, over_boundary_record)

        pos_bak = pos.copy() #pre_pos = xdatcar[n_step-1].get_positions()
        # 校准位置
        pos = calibration_pos(pos,xdatcar_0.cell,over_boundary_record)
        #over_boundary_record = np.zeros((n_atoms_total,3)) #初始化
        
        msd_list_append = []
        msd_list_append.append((n_step+1)*optim) #输出当前时间
        msd_list_append.append(calc_msd(pos,xdatcar_0.get_positions(),3)) #输出总msd
        msd_list_append.extend(calc_msd(pos,xdatcar_0.get_positions(),1)) #输出pmsd
        msd_list_append.extend(calc_msd(pos,xdatcar_0.get_positions(),3,list(atoms_info_dict.values()))) #输出pmsd
        msd_list  = np.vstack([msd_list,msd_list_append])
        
        #if n_step>=10: break
    np.savetxt("msd.out",msd_list,fmt="%s")
    plot()  #绘图
    
    time_cost = time.perf_counter () - start_time
    print("---- Time used: %s ----" % time.strftime("%H:%M:%S", time.gmtime(time_cost)))


