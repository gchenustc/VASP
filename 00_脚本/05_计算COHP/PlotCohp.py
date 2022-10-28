import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# ******** 配置 ********
plt.rcParams['font.sans-serif'] = ['Arial'] # 字体
cohp_x = "auto" # COHP x 轴范围, 可选"auto"
cohp_y= [-27, 10]# COHP y 轴范围，可选"auto"
cohp_sep_x = 1 # x轴刻度间隔
cohp_sep_y = 4 # y轴刻度间隔
cohp_0_plus_color = "blue"# -cohp大于0的颜色，可以用 #123456 表示
cohp_0_plus_alhpa = 0.2 # -cohp大于0的透明度
cohp_0_minus_color = 'red'# -cohp小于0的颜色
cohp_0_minus_alhpa = 0.2 # -cohp小于0的透明度
icohp_color = 'green'# icohp显示的颜色
# ******** END 配置 ********

def read_COHP(fn):
    """
    return.shape[0]指的是cohp的能量范围, return.shape[1]是计算的cohp的数目
    """
    raw = open(fn).readlines()
    
    raw = [l for l in raw if 'No' not in l][3:]
    raw = [[eval(i) for i in l.split()] for l in raw]  # 将字符串转为数字
    return np.array(raw)

data_cohp = read_COHP('./COHPCAR.lobster')
labels_cohp = [l[:-1] for l in open('./labels').readlines()]  # l[:-1]对的目的是去除 \n
len_cal = len(labels_cohp)
icohp_ef = [eval(l.split()[-1])
            for l in open('./ICOHPLIST.lobster').readlines()[1:]]  # [-11.53014,-11.44923,-11.40647,...]

# for i in range(len_cal):
for i in range(1):
    fig, ax1 = plt.subplots(figsize=([2.4, 4.8]),dpi=100)
    # ************ COHP ************
    ax1.plot(-data_cohp[:, i*2+3], data_cohp[:, 0], color='k', label=labels_cohp[i]) # 因为画的是-COHP，x轴加负号
        # cohp fill > 0
    ax1.fill_betweenx(y=data_cohp[:, 0], x1=-data_cohp[:, i*2+3], x2=0, where=-data_cohp[:, i*2+3] >= 0, facecolor=cohp_0_plus_color, alpha=cohp_0_plus_alhpa)
        # cohp fill < 0
    ax1.fill_betweenx(y=data_cohp[:, 0], x1=-data_cohp[:, i*2+3], x2=0, where=-data_cohp[:, i*2+3] <= 0, facecolor=cohp_0_minus_color, alpha=cohp_0_minus_alhpa)
        # cohp 的显示范围
    if isinstance(cohp_x, str) and cohp_x.lower()=="auto":
        pass
    else:
        ax1.set_xlim(cohp_x)  # 根据lobster的设置更改, lobster算了某个范围的能量后才有效果
    if isinstance(cohp_y, str) and cohp_y.lower()=="auto":
        pass
    else:
        ax1.set_ylim(cohp_y)  # 这一行随意
        # x轴y轴label
    ax1.set_xlabel('-COHP (eV)', color='k', fontsize='large')
    ax1.set_ylabel('$E-E_\mathrm{F}$ (eV)', fontsize='large')
        # tick
            # x 频率
    ax1.xaxis.set_major_locator(MultipleLocator(cohp_sep_x))  # x轴每隔多少输出一个刻度
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2)) # x轴每两个大刻度之间有几个小刻度，1代表无
            # y 频率
    ax1.yaxis.set_major_locator(MultipleLocator(cohp_sep_y))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
            # tick 配置
    ax1.tick_params(which='major', labelsize=9, labelcolor='k', direction='in')
    ax1.tick_params(which='minor', labelsize=9, labelcolor='k', direction='in')
    # ************ END COHP ************
    
    # ************ ICOHP ************
    ax2 = ax1.twiny()
    ax2.plot(-data_cohp[:, i*2+4], data_cohp[:, 0], color=icohp_color)
        # icohp 的显示范围, 推荐注释掉，线条会自动调整
    #ax2.set_ylim([-30, 10])
    #ax2.set_xlim([-0.01, 18])
    ax2.set_xlabel('-ICOHP (eV)', color=icohp_color, fontsize='large')
        # tick 和 tick label 的位置
    ax2.xaxis.tick_top()  # tick_bottom()
        # axis label的位置
    ax2.xaxis.set_label_position('top')
        # tick
            # x 频率
    ax2.xaxis.set_major_locator(MultipleLocator(3))  # x轴每隔多少输出一个刻度
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2)) # x轴每两个大刻度之间有几个小刻度，1代表无
            # y 频率 
    #ax2.yaxis.set_major_locator(MultipleLocator(5))
    #ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
            # tick 配置
    ax2.tick_params(which='major', labelsize=9, labelcolor='k', direction='in')
    ax2.tick_params(which='minor', labelsize=9, labelcolor='k', direction='in')
    
    # ************ END ICOHP ************

    # ************ markers ************
    ax1.axvline(x=0, color='k', linestyle=':', alpha=0.5)
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.5)
#     ax2.annotate(labels_cohp[i], xy=(17.5, 9.5), ha='right', va='top', bbox=dict(boxstyle='round', fc='w', alpha=0.3))
    ax1.legend(loc='upper right', borderaxespad=0.5, borderpad=0.4, fontsize=9, frameon=True, fancybox=True, framealpha=0.6, handlelength=0, handleheight=0, handletextpad=0)
    ax2.annotate(f'{-icohp_ef[i]:.3f}', xy=(17.5, -0.08),
                 ha='right', va='top', color='blue')
    fig.savefig(f'cohp-{i+1}.png', dpi=500,
                bbox_inches="tight", transparent=True)
    #plt.close()
    # ************ end markers ************