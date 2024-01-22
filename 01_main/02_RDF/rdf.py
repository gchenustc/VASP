"""
Desc: 
    1. Calculating the rdf for molecular dynamics (vasp). the first and the second params are the frame indices for the calculated range. the third param is the max distance for the rdf (rmax). this value less than the half of the shortest axis' length. the fourth param is the number of division for the rmax. the fifth param is the appointed element.  
    2. eg: if calculate the N atom's rdf from no.1 to no.3000 frames, 3.0 is the half of the shortest axis for cell and length's division equals 70. Type on the command line: python script.py 1 3000 3.0 70 N.
Author: 
    gchen
Time: 
    20240118
"""

import sys
import numpy as np
from ase.io.vasp import read_vasp_xdatcar
from ase.geometry.analysis import Analysis
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import argrelextrema
from scipy.interpolate import interp1d
mpl.rcParams['font.family'] = 'Arial'


####### setting ######
file_name = "rdf.csv"
fig_name = "rdf.png"
smooth = False  # if use the smooth.
# the higher of the ws, the smoother of the line. this param is available after set the true for the smooth param.
ws = 6
####### setting ######

n1 = int(sys.argv[1])-1  # the start frame of XDATCAR
n2 = int(sys.argv[2])-1  # the end frame of XDATCAR
# the max distance for the rdf, less then the half of the shortest axis
rmax = float(sys.argv[3])
# the number of division for the rmax
nbins = int(sys.argv[4])
element = sys.argv[5]

# poscar = read("POSCAR")
xdatcar = read_vasp_xdatcar("XDATCAR", index=slice(0, None, 1))

ana = Analysis(xdatcar)

count = n2-n1+1  # the number of frames.
rdf = 0
for i in range(n1, n2+1):
    tmp = ana.get_rdf(rmax=rmax, nbins=nbins, imageIdx=i, elements=element, return_dists=True)[
        0]  # tmp[0] if rdf value, tmp[1] if the distance
    # print(tmp[1])
    rdf += tmp[0]
rdf /= count


def plotrdf(x, y, name, smooth=True, window_size=5):
    # 创建画布和子图
    fig, ax = plt.subplots(figsize=(6, 4), dpi=100)

    if smooth:
        # 使用移动平均滤波器平滑数据
        y = np.convolve(y, np.ones(window_size)/window_size, mode='same')
        # 使用线性插值方法平滑数据
        # interp_func = interp1d(x, y, kind='linear')
        # x = np.linspace(x.min(), x.max(), 3000)
        # y = interp_func(x)

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
    # ax.legend()

    # 保存图片
    plt.savefig(name, dpi=300)

    # 打印峰谷
    maxima = argrelextrema(y, np.greater, order=1)[0]  # 谷峰的索引
    for i in maxima:
        print("峰谷:", "x:", x[i], "y:", y[i])


def write_csv(x, y, name):
    with open(name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["dist(A)", "rdf"])
        writer.writerows(zip(x, y))


def main():
    plotrdf(tmp[1], rdf, fig_name, smooth=smooth, window_size=ws)
    write_csv(tmp[1], rdf, file_name)


if __name__ == "__main__":
    main()