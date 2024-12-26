import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


##############################

# 根据MD模拟确定一些参数

dt = 0.001  # 积分步长，单位为 ps
T = 300  # 单位为 K
a = 30.4219 # 晶格常数，单位为 A

font_size = 12 # 设置字体大小

##############################

# 常量
kB = 8.617e-5  # 玻尔兹曼常数，单位为 eV/K
V = (4 * a) ** 3  # 体积

##############################

# 读取文件
file = 'hac.txt'
# 读取前两列以获取Tstart和Nc (Nc: 一次模拟的关联数据量)
with open(file, 'r') as f:
    for _ in range(3):  # 跳过前三行
        next(f)
    Tstart, Nc = map(int, f.readline().split())

print(f"Np*step (time of one turn (fs)): {Tstart}, Nc (samples number of one turn): {Nc}")

# 读取后续数据 (Ns = 取样间隔)
data = pd.read_csv(file, header=None, sep='\s+', skiprows=Nc+4, names=['index', 'TimeDelta', 'Ncount', 'c_1', 'c_2', 'c_3'])

# print(data)

Ns = data['TimeDelta'].iloc[2] - data['TimeDelta'].iloc[1]  # 取样间隔
print(f"Ns (the sample interval): {Ns}")

# 计算缩放系数
scale = 1.6e+3 / ((kB * T * T * V) * (Ns * dt))  # 转换因子

# 创建hac数组并移除每Nc+1行的第一个元素
hac = data.values
# print(hac)
nall = len(hac) // (Nc + 1)

# 删除每Nc+1行的第一个元素
hac = np.delete(hac, np.arange(0, nall * (Nc + 1), Nc + 1), axis=0)

# 清除不必要的变量
del data

# M: 自动确定模拟次数
M = len(hac) // Nc
print(f"M: (number of turn): {M}")

# 关联时间
t = np.arange(1, Nc + 1) * dt * Ns

# 计算三维系统的热流自关联和跑动热导率
# 平均热流
hac = np.mean(hac[:, 3:6], axis=1).reshape(M, Nc)
hac = np.transpose(hac)
print(f"hac: {hac}, hac.shape: {hac.shape}")

rtc = scale * np.cumsum(hac, axis=0) - 0.5 * scale * hac  # 使用cumtrapz的近似方法
print(f"rtc: {rtc}, rtc.shape: {rtc.shape}")

# 系综平均（对独立模拟的平均）
hac_ave = np.mean(hac, axis=1)
rtc_ave = np.mean(rtc, axis=1)
print(f"hac_ave: {hac_ave}")
print(f"rtc_ave: {rtc_ave}")

# 绘制归一化的热流自关联
plt.figure()
for n in range(M):
    plt.plot(t, hac[:, n] / hac[0, n], color='0.5')
plt.plot(t, hac_ave / hac_ave[0], linewidth=3)
plt.xlabel('Correlation Time (ps)', fontsize=font_size)
plt.ylabel('HCACF (Normalized)', fontsize=font_size)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.ylim(-0.7, 1.3)
plt.savefig("Heat_flo_autocorrelation.png")
plt.show()

# 绘制跑动热导率
plt.figure()
for n in range(M):
    plt.plot(t, rtc[:, n], color='0.5')
plt.plot(t, rtc_ave, linewidth=3)
plt.xlabel('Correlation Time (ps)', fontsize=font_size)
plt.ylabel(r'$\kappa$ (W/mK)', fontsize=font_size)
plt.ylim([100, 600])
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("dynamic_thermal_conductivity.png")
plt.show()

# 报道结果
kappa_converged = np.mean(rtc[Nc//2:], axis=0)  # 具体问题具体分析
kappa_average = np.mean(kappa_converged)
kappa_error = np.std(kappa_converged) / np.sqrt(M)

print(f'Converged kappa: {kappa_average}')
print(f'Error in kappa: {kappa_error}')
