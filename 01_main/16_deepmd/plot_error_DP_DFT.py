import matplotlib.pyplot as plt
import numpy as np
"""
Author: gchen
Date: 20241113
Description: 
    该脚本用于从指定文件中(detile.f.out, 该文件由dp test命令产生)读取vasp和dp的原子受力信息，并使用 matplotlib 绘制散点图和 y=x 的直线。
    同时，从日志文件中(out.log)读取 Force MAE 和 Force RMSE 的值，并在图例中显示这些指标。
    detile.f.out和out.log需要dp test生成，比如: dp test -m demo.graph -s data -d detile -l out.log

Usage:
    1. 将数据文件 `detile.f.out` 和日志文件 `out.log` 放在脚本同一目录下。
    2. 运行脚本：`python script_name.py`
    3. 图表将保存为 `f.png` 文件。

Parameters:
    - force_file_path: 原子受力数据文件路径，默认为 'detile.f.out'
    - log_file_path: 日志文件路径，默认为 'out.log'
"""

def read_data(file_path):
    """
    从给定的文件路径(detile.f.out)读取数据，返回 x 和 y 数据列表。
    
    参数:
    - file_path: 包含数据的文件路径
    
    返回:
    - x_data: x 轴数据列表
    - y_data: y 轴数据列表
    """
    
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]


    x_data = []
    y_data = []

    for line in lines:
        if "#" in line:
            continue
        columns = line.strip().split()  # 假设数据是以空格分隔的
        x_data.append(float(columns[0]))
        x_data.append(float(columns[1]))
        x_data.append(float(columns[2]))
        y_data.append(float(columns[3]))
        y_data.append(float(columns[4]))
        y_data.append(float(columns[5]))

    return x_data, y_data

def read_metrics(log_file):
    """
    从dp test的日志文件中读取 Force MAE 和 Force RMSE 的值。dp test需要给定指定参数（dp test -m demo.graph -s data -l out.log）生成日志文件，其中out.log为日志文件。
    
    参数:
    - log_file: 包含日志信息的文件路径
    
    返回:
    - mae: Force MAE 的值字符串
    - rmse: Force RMSE 的值字符串
    """
    mae = ""
    rmse = ""

    with open(log_file, 'r') as file:
        for line in file:
            if 'Force  MAE' in line:
                mae += line.split('  : ')[1].strip() + "; "
            elif 'Force  RMSE' in line:
                rmse += line.split('  : ')[1].strip() + "; "

    return mae, rmse

def plot_data(x_data, y_data, mae, rmse):
    """
    使用 matplotlib 绘制点图和 y=x 的直线。
    
    参数:
    - x_data: x 轴数据列表: vasp计算的力
    - y_data: y 轴数据列表：dp预测的力
    - mae: Force MAE 的值字符串
    - rmse: Force RMSE 的值字符串
    """
    
    # 创建一个新的图形
    plt.figure(figsize=(8, 8))

    # 绘制点图
    plt.scatter(x_data, y_data, color='blue', label=f'MAE: {mae}\nRMSE: {rmse}')

    # 绘制 y=x 的直线
    #lims = [np.min([plt.xlim(), plt.ylim()]), np.max([plt.xlim(), plt.ylim()])]
    plt.plot([-6,6], [-6,6], color='red', linestyle='--', linewidth=3)

    # 设置图表属性
    plt.xlabel('vasp force (eV/A)', fontsize=16)
    plt.ylabel('pred force (eV/A)', fontsize=16)
    plt.title('')
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    #plt.axis('equal')  # 保持纵横比一致

    # 显示图表
    plt.savefig('f.png')

# 文件路径
force_file_path = 'detile.f.out'
log_file_path = "out.log"

# 读取数据
x_data, y_data = read_data(force_file_path)

# 读取力的平均误差和均方误差
mae, rmse = read_metrics(log_file_path)

# 绘制数据
plot_data(x_data, y_data, mae, rmse)
