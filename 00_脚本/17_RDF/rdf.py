from ase.io.vasp import read_vasp_xdatcar
from ase.io import read
from ase.geometry.analysis import Analysis
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.interpolate import interp1d
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Arial'
"""
绘制RDF
"""

# 设定数据 ######
n = 1   #绘制第多少帧的rdf
rmax=3.3   # rdf最大距离，要小于最短轴的一半
#interval=0.0002   # 距离检测间隔
fig_name = "rdf.png"   # 保存图片的名称
ws=1    # 线条的平滑程度，越高越平滑，但也越不准确
nbins=60   # 密度，不能太大，否则在真空处RDF会无限大，会导致线条不连续
##################
n=n-1

# 导入数据
#poscar = read("POSCAR")
xdatcar = read_vasp_xdatcar("XDATCAR",index=slice(0,None,1))

ana=Analysis(xdatcar)

def plotrdf(x,y,name,smooth=True, window_size=5):
    # 创建画布和子图
    fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
    
    if smooth:
        # 使用移动平均滤波器平滑数据
        y = np.convolve(y, np.ones(window_size)/window_size, mode='same')
        # 使用线性插值方法平滑数据
        #interp_func = interp1d(x, y, kind='linear')
        #x = np.linspace(x.min(), x.max(), 3000)
        #y = interp_func(x)

    # 绘制线形图
    ax.plot(x, y, linestyle='-', color='b', label='Line')
    

    # 设置标题和坐标轴标签
    ax.set_title('RDF')
    ax.set_xlabel('dist(A)')
    ax.set_ylabel('rdf')
    
    # 设置字体大小
    ax.title.set_fontsize(18)
    ax.xaxis.label.set_fontsize(15)
    ax.yaxis.label.set_fontsize(15)

    # 设置图例
    #ax.legend()

    # 保存图片
    plt.savefig(name, dpi=300)
    
    #
    maxima = argrelextrema(y, np.greater, order=1)[0] # 谷峰的索引
    
    for i in maxima:
        print("峰谷",x[i],y[i])
        
ret = ana.get_rdf(rmax=rmax, nbins=nbins, imageIdx=n, elements=None, return_dists=True)[0]
plotrdf(ret[1],ret[0], fig_name, smooth=False, window_size=ws)


with open('rdf.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["dist(A)","rdf"])
    writer.writerows(zip(ret[1],ret[0]))
