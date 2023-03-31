import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['font.sans-serif']="Arial"

 
#输入因变量
x = np.array(["0", "20", "40", "60", "80" ,"100"])

y1 = np.array([1.4067, 1.3870, 1.3721, 1.3598, 1.3495, 1.3406])
y2 = np.array([1.356, 1.333, 1.319, 1.307, 1.297, 1.288])
y3 = np.array([1.32784, 1.31369, 1.30601, 1.29755, 1.29027, 1.28385])
#y4 = np.array([0.879261, 0.770276, 0.893485, 0.892955, 0.892227, 0.890149])
ys=[y1,y2,y3]
markers=["s","v","o","D"]
#linestyles=["-","--","-.",":"]
linestyles=["-","-","-","-"]
labels=["cg-N","N$_{3}^{-}$","N$_{5}^{-}$"]

 
fig,ax=plt.subplots(figsize=(6.4,4.8), dpi=300)
#设置自变量的范围和个数
#画图
for index,i in enumerate(ys):
    ax.plot(x,i, label=labels[index], linestyle=linestyles[index], marker=markers[index], markersize='6', linewidth=1.5)
#设置坐标轴
ax.set_xlim(0, 5)
ax.set_ylim(1.23, 1.47)
ax.set_xlabel('压力 (GPa)', fontsize=13)
ax.set_ylabel('键长 (Å)', fontsize=13)
#设置刻度
ax.minorticks_on()
ax.tick_params(axis='both', which="major", width=1, labelsize=11, direction="in")
ax.tick_params(axis='both', which="minor", width=0.6, labelsize=11, direction="in")
# 设置边框宽度
ax.spines['left'].set_linewidth(1.2)
ax.spines['right'].set_linewidth(1.2)
ax.spines['top'].set_linewidth(1.2)
ax.spines['bottom'].set_linewidth(1.2)
# 绘制参考线和注释
ax.axhline(1.10, c='black', ls='--', lw=1)
#ax.text(x=5.30,y=1.07,s="三键", size=10, alpha=0.8)
ax.axhline(1.25, c='black', ls='--', lw=1)
#ax.text(x=5.30,y=1.22,s="双键", size=10, alpha=0.8)
ax.axhline(1.45, c='black', ls='--', lw=1)
#ax.text(x=5.30,y=1.42,s="单键", size=10, alpha=0.8)
# ax.axvline(1)
#显示网格
#ax.grid(True, linestyle='-.')
#ax.yaxis.grid(True, linestyle='-.')
#添加图例
legend = ax.legend(frameon=False,loc='best')
 
# plt.show()
fig.savefig('out.png')
