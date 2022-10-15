import numpy as np
import pandas as pd
import os
import re
import math
import seaborn as sns
import matplotlib.pyplot as plt

# 获得某一步晶胞常数,numpy.array形式返回
def get_onestep_constant(onestep_data):
    head_info = re.split(r"Direct configuration= +\d+",onestep_data)[0]
    xyz_coord = np.array(re.findall(r"-?\d+\.[\d]{1,}",head_info)).reshape(3,3).astype(np.float64) #一共九个坐标
    return xyz_coord

# 给定一个原子的分数坐标和该胞的晶格常数，得到笛卡尔坐标
def get_cartesol_coord(fraction_coord,cell_constant):
    cartesol_coord=[]
    for index in range(3):
        ss = fraction_coord[0]*cell_constant[0][index]+fraction_coord[1]*cell_constant[1][index]+fraction_coord[2]*cell_constant[2][index]
        cartesol_coord.append(ss)
    return cartesol_coord


# 获得某一步的原子的所有坐标
def get_onestep_atoms_coord(onestep_data):
    atoms_coord_str = re.split(r"Direct configuration= +\d+\n",onestep_data)[1]
    atoms_coord = re.split(r"\n",atoms_coord_str)[:-1]
        # 去除每一行原子坐标左右的空格,并且将以字符串打印出来的具体坐标切分成以x,y,z轴分割的列表列表
    for atom_index in range(len(atoms_coord)):
        atoms_coord.insert(atom_index, re.split(r"\s+",atoms_coord.pop(atom_index).strip()))
        # 将str转为int格式
    atoms_coord = np.array(atoms_coord).astype(np.float64)
        # 转为笛卡尔坐标
    atoms_coord_cartesol=[]
    cell_constant = get_onestep_constant(onestep_data)
    for atom_coord in atoms_coord:
        atoms_coord_cartesol.append(get_cartesol_coord(atom_coord,cell_constant))
    return atoms_coord_cartesol

# 获得某一步的指定原子的坐标,designated参数为指定的编号的原子
def get_onestep_atoms_coord_designated(onestep_data,*designated):
    atoms_coord = get_onestep_atoms_coord(onestep_data)
    stoms_coord_designated = []
        #如果传入的参数是列表的话，则使用列表
    if isinstance(designated[0],list):designated=designated[0]
    for index in designated:
        stoms_coord_designated.append(atoms_coord[index-1])
    return stoms_coord_designated

# 获得x轴步数和y轴距离
def get_step_dist(designated_atoms, steps_to_start_plot, file_name, center_atoms):
    with open(file_name,"r") as f:
        data = f.read()
    data_list = re.split(r"^[A-CE-Z].+\n",data,flags=re.M)[1:]

    # step为一步的XDATCAR信息
    y_axis_distance = []
    for step in range(len(data_list)):
        atoms_coord = get_onestep_atoms_coord_designated(data_list[step],*designated_atoms)
#         center_coord = get_cartesol_coord([0.5,0.5,0.5],get_onestep_constant(data_list[step]))
        
        if center_atoms:
            center_coord_np = np.array(get_onestep_atoms_coord_designated(data_list[step],center_atoms))
        else:  # 中心原子的坐标为中间一层
            center = len(get_onestep_atoms_coord(data_list[step]))//2
            center_coord_np = np.array(get_onestep_atoms_coord_designated(data_list[step],[i for i in range(center-4,center+5)]))
        center_coord_pd = pd.DataFrame(center_coord_np)
        center_coord = [center_coord_pd[0].mean(),center_coord_pd[1].mean(),center_coord_pd[2].mean()]
        # 指定原子距离中心点的距离的列表
        distance=[]
        for atom_coord in atoms_coord:
            distance.append(math.sqrt(math.pow(atom_coord[0]-center_coord[0],2)+\
                                    math.pow(atom_coord[1]-center_coord[1],2)+\
                                    math.pow(atom_coord[2]-center_coord[2],2)))
        distance_mean = np.array(distance).mean() # 平均距离
        y_axis_distance.append(distance_mean)
    y_axis_distance = y_axis_distance[steps_to_start_plot-1:]
    dataframe = pd.DataFrame({"step":range(steps_to_start_plot,len(y_axis_distance)+steps_to_start_plot),"dist":y_axis_distance})
    sns.lineplot(x="step",y="dist",data=dataframe,ax=axes)

    #显示所有列
    pd.set_option('display.max_columns', None)
    #显示所有行
    pd.set_option('display.max_rows',None)
    print(dataframe)


if __name__ == "__main__":
    file_name = input("请输入文件名")
    line_num = int(input("在一个图上绘制几条线？"))
    steps_to_start_plot = int(input("请输入从第几步开始绘图:"))
    if input("需要手动确定中心吗,需要请输入任意字符"):
        center_atoms = list(np.array(input("请输入中心原子，相邻两个原子序号用空格切分").split(" ")).astype(np.int8))
    else: center_atoms=""
    fig,axes = plt.subplots(1,1,figsize=(15,10))
    for i in range(line_num):
        print("绘制第"+str(i+1)+"条：")
        designated_atoms = np.array(input("请输入原子序号，相邻两个原子序号用空格切分:").split(" ")).astype(np.int8)
        get_step_dist(designated_atoms, steps_to_start_plot, file_name, center_atoms)
    fig.savefig("fig.png",dpi=400)
