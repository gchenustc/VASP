import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.io import read
from ase.units import Bohr, Angstrom
"""
Author: gchen
Time: 20250120
Description: 这段代码是一个用于分析分子动力学模拟中均方位移（MSD）的工具，它包括了从VASP软件生成的XDATCAR文件读取轨迹数据、计算每个原子类型的MSD、绘制MSD图以及导出结果为CSV文件的功能。具体来说，该脚本实现了以下功能：
读取轨迹：通过read_vasp_trajectory函数读取XDATCAR文件中的所有结构快照，并提取原子位置和化学符号。
计算MSD：使用compute_msd函数计算整个系统的MSD以及每个元素在不同方向上的MSD。
绘制图表：通过plot_msd函数生成MSD随时间变化的图表，并保存为图片文件。
导出数据：利用export_data函数将MSD的数据导出为CSV格式，以便进一步分析或存档。
Usage: 确保你有一个由VASP生成的XDATCAR文件，里面包含有分子动力学模拟的时间序列快照。同时，你需要安装Python环境及依赖库（如numpy, pandas, matplotlib, 和 ase），这些库可以通过pip进行安装：pip install numpy pandas matplotlib ase
1. 将此代码保存为一个Python脚本文件，例如msd_analysis.py。
2. 打开命令行界面，导航到存放msd_analysis.py和XDATCAR文件的目录。
3. 运行脚本：python msd_analysis.py
4. 当提示输入时键入模拟步长的时间间隔（timestep），以飞秒（fs）为单位。

运行后，该脚本会执行以下操作：

计算并输出一个名为msd_output.csv的CSV文件，其中包含了所有计算出来的MSD数据。
生成三个PNG图像文件：
msd.png：显示系统整体的MSD随时间的变化情况。
msd_elements.png：显示各个元素总的MSD随时间的变化情况。
msd_p_elements.png：显示各个元素在x, y, z方向上MSD随时间的变化情况。
"""

def read_vasp_trajectory(filename='XDATCAR'):
    """Read VASP trajectory file and return positions and symbols."""
    traj = read(filename, index=':')
    symbols = traj[0].get_chemical_symbols()
    positions = [atoms.get_positions() for atoms in traj]
    unique_symbols = list(set(symbols))
    return np.array(positions), symbols, unique_symbols

def compute_msd(positions, symbols, timestep):
    """Compute the mean squared displacement (MSD) for each atom type and direction."""
    unique_symbols = list(set(symbols))
    n_steps, n_atoms, _ = positions.shape
    msd_data = {symbol: {'total': [], 'x': [], 'y': [], 'z': []} for symbol in unique_symbols}
    total_msd = {'total': [], 'x': [], 'y': [], 'z': []}
    time_data = []

    for step in range(1, n_steps):
        time_data.append((step+1)*timestep)
        disp = positions[step] - positions[0]
        msd_total = np.sum(disp ** 2, axis=1)
        msd_x = disp[:, 0] ** 2
        msd_y = disp[:, 1] ** 2
        msd_z = disp[:, 2] ** 2

        for symbol in unique_symbols:
            mask = np.array([s == symbol for s in symbols])
            msd_data[symbol]['total'].append(np.mean(msd_total[mask]))
            msd_data[symbol]['x'].append(np.mean(msd_x[mask]))
            msd_data[symbol]['y'].append(np.mean(msd_y[mask]))
            msd_data[symbol]['z'].append(np.mean(msd_z[mask]))

        total_msd['total'].append(np.mean(msd_total))
        total_msd['x'].append(np.mean(msd_x))
        total_msd['y'].append(np.mean(msd_y))
        total_msd['z'].append(np.mean(msd_z))

    # 创建一个空的 DataFrame 来存放所有的数据
    data = {'time': time_data}
    for key in ['total', 'x', 'y', 'z']:
        data[key] = total_msd[key]

    # 添加每个 symbol 的数据
    for symbol in unique_symbols:
        for key in ['total', 'x', 'y', 'z']:
            column_name = f'{symbol}_{key}'
            data[column_name] = msd_data[symbol][key]

    return data

def plot_msd(data, data_labels, output_file='msd.png'):
    df = pd.DataFrame(data)

    """Plot MSD data for each element and direction."""
    fig, ax = plt.subplots(figsize=(10, 6))
    #steps = np.arange(1, len(total_msd['total']) + 1)

    # 绘制所有元素综合的MSD
    for key in data_labels:
        ax.plot(df['time'], df[key], label=f'{key}')

    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('MSD (Å²)')
    ax.set_title('Mean Squared Displacement (MSD)')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlim(0,data['time'][-1])
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()
  

def export_data(data):
    # 将字典转换成 DataFrame
    df = pd.DataFrame(data)

    # 导出为 CSV 文件
    df.to_csv('msd_output.csv', index=False)


if __name__ == "__main__":
    # User-defined timestep in femtoseconds (fs)
    timestep = float(input("Enter the timestep (in fs): "))
    
    # Read trajectory data
    positions, symbols, unique_symbols = read_vasp_trajectory('XDATCAR')

    # Compute MSD
    data = compute_msd(positions, symbols, timestep)
    
    export_data(data)

    # Plot MSD
    plot_msd(data, ['total', 'x', 'y', 'z'], 'msd.png')
    
    element_labels = []
    for key in unique_symbols:
        element_labels.append(key+'_total')
    plot_msd(data, element_labels, "msd_elements.png")

    pelement_labels = []
    for key in unique_symbols:
        pelement_labels.append(key+'_x')
        pelement_labels.append(key+'_y')
        pelement_labels.append(key+'_z')
    plot_msd(data, pelement_labels, "msd_p_elements.png")
    
   
 
    