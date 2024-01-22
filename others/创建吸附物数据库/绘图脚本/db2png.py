import ase
from ase.io import read
from ase.db import connect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib import rcParams 
rcParams['font.family']='Arial'


def selectAtomsRows(mydb, temp_info_list=None):
    """
    temp_info 是模板结构的信息，是个列表 --> [[n_atoms, n_vac, energy],[...],...]
    """
    energy_H = -3.38794825 # 用DFT计算的氢气的能量/2
    
    selected_rows_list = {}  #{[96,8,-650.074]: [row1, row2, ...], [...]:[row1,row2]}
    for i in temp_info_list:
        n_atoms, n_vac, energy = i
        selected_rows_list[tuple(i)] = list(mydb.select(f"N={n_atoms}, H>0, relaxed=True, relaxed_converg=True"))
    
    select_rows_info = []  # ["id": 1, "energy":-960, "natoms":98, "n_N":90, "n_H"":8, "n_vac":10, adsorp_ratio":0.8, "adsorp_e":-0.02]
    for key,rows in selected_rows_list.items():
        for row in rows:
            dict_ = {}
            dict_["id"] = row.id
            dict_["energy"] = row.energy
            dict_["natoms"] = row.natoms
            dict_["n_N"] = key[0]
            n_adsorps = np.sum(row.numbers==1) # 1是吸附物H的id
            dict_["n_H"] = n_adsorps
            adsorp_e = (row.energy - key[2] - energy_H * n_adsorps)/n_adsorps
            dict_["n_vac"] = key[1]
            #dict_["adsorp_ratio"] = float(str(n_adsorps/key[1])[:5])
            dict_["adsorp_ratio"] = "{:.2f}".format(n_adsorps/key[1])
            dict_["adsorp_e"] = adsorp_e
            #print(dict_)
            select_rows_info.append(dict_)
            
    return select_rows_info


def plot(selected_rows_info):
    x_axis = [i["adsorp_ratio"] for i in selected_rows_info]
    y_axis = [i["adsorp_e"] for i in selected_rows_info]
    df = pd.DataFrame({"adsorp_ratio":x_axis, "adsorp_energy":y_axis})
    print(df)
    df_counts = df.groupby(['adsorp_ratio', 'adsorp_energy']).size().reset_index(name='counts')
    # Draw Stripplot
    fig, ax = plt.subplots(figsize=(6,5), dpi= 150)
    sns.stripplot(x=df_counts.adsorp_ratio, y=df_counts.adsorp_energy, sizes=df_counts.counts*200, ax=ax, edgecolor="gray", linewidth=1, jitter=0.2) # jitter: 抖动量

    # Decorations
    plt.title('100a Adsorption Energy', fontsize=22)
    #plt.show()
    #plt.axhline(0, color ="black", linestyle ="--", lw=2)
    
    plt.savefig(f"{db_name}.png")

db_name = "111a-vac-H.db"
mydb = connect(db_name)
row_infos = selectAtomsRows(mydb, temp_info_list=[[96,12,-641.689]]) # 二级列表：氮原子数，H空位数，纯氮的能量；一级列表：数据库中有多种扩胞方式的结构，则填写多个：比如[[96,12,-641.689],[144,18,-888]]
plot(row_infos)