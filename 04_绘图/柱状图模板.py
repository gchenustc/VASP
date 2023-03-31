import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import warnings
warnings.filterwarnings("ignore")

"""
from matplotlib.font_manager import FontManager  
fonts = set(f.name for f in FontManager().ttflist)  
print ('可用字体:')  
for f in sorted(fonts):  
    print('\t' + f)
"""

plt.rcParams['font.family'] = ['Arial'] # SimHei
fig,ax = plt.subplots(1,1,figsize=(6.4,4.8),dpi=300)

labels = ["cg-N",
"LP-N",
"HLP-N",
"BP-N",
"LiN$_5$",
"NaN$_5$",
"CsN$_5$",
"K$_2$N$_6$",
"K$_9$N$_{56}$"
]

x_kind1 = [50, 188, 280, 0, 9.9, 20, 9.1, 40, 0]
x_kind2 = [110, 150, 244, 146, 45, 50, 60, 45, 46]
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

label_font = {
    #'weight':'bold',
    'size':14,
    'family':'Arial'
}

rects1 = ax.bar(x - width/2, x_kind1, width, label='理论压力',ec='k',color="#709bff",lw=.8
               )
rects2 = ax.bar(x + width/2 + .05, x_kind2, width, label='实验压力',ec='k',color='#ff5757',
                lw=.8)
# rects2 = ax.bar(x + width/2 + .05, means_2016, width, label='2016',ec='k',color='k',
#                 lw=.8,alpha=.8)
ax.tick_params(which='major',direction='in',length=5,width=1.5,labelsize=11,bottom=False)
ax.tick_params(axis='x',labelsize=11,bottom=False,labelrotation=15)
ax.set_xticks(x)
ax.set_ylim(ymin = 0,ymax = 300)
ax.set_yticks(np.arange(0,300,30))

ax.set_ylabel('压力 (GPa)',fontdict=label_font)
ax.set_xticklabels(labels,fontdict=label_font)
#ax.legend(markerscale=10,fontsize=12,prop=legend_font)
ax.legend(markerscale=10,fontsize=12,frameon =False)

# Add some text for labels, title and custom x-axis tick labels, etc.
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{:.0f}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
autolabel(rects2)
fig.tight_layout()
plt.savefig(r'barplot05.png',bbox_inches = 'tight')