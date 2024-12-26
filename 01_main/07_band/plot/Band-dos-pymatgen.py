# plot_vasp_band_structure.py
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter, DosPlotter, BSPlotter
import matplotlib.pyplot as plt
#from matplotlib import font_manager 
#font_path = "/public/home/c1337375425/storage/fonts"
#font_files = font_manager.findSystemFonts(fontpaths=font_path)
#for file in font_files:
#    font_manager.fontManager.addfont(file)
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# 变量初始化
bs_vasprun = Vasprun("vasprun.xml", parse_projected_eigen=True)	
bs_data = bs_vasprun.get_band_structure(line_mode=True)  # 杂化泛函需要同时打开：line_mode=True, force_hybrid_mode=True，需要注意的是杂化泛函的能带计算也要读入k点信息，所以保证在KPOINTS文件里的高对称k点的后面添加标记，这样get_band_structure函数可以自动的识别这些高对称点，否则，需要手动载入包含k点信息的字典。
dos_data = bs_vasprun.complete_dos

"""
BSDOSPlotter的参数解释：cb_energy_range：导带从有带的地方开始，能量向上画到结束的区间。vb_energy_range：同理
"""
vb_energy_range = 8     # 从开始出现价带算起最大的y轴范围
cb_energy_range = 8     # 从开始出现导带算起最大的y轴范围

def plot_band():
    bands_plot = BSPlotter(bs = bs_data)
    bands_plot.get_plot(ylim=[-vb_energy_range, cb_energy_range]) 
    # 如果需要保存图片
    # 支持格式 png, tif, jpg....，只要改后缀名即可
    plt.savefig("band.png")

def plot_dos():
    dos_plot = DosPlotter()
    dos_plot.add_dos("Total DOS", dos_data)
    dos_plot.get_plot([-vb_energy_range, cb_energy_range], [0, 12])
    plt.savefig("dos.png")

def plot_bands_dos():
    bs_dos_plot = BSDOSPlotter(bs_projection='elements', dos_projection='elements', vb_energy_range=vb_energy_range, cb_energy_range=cb_energy_range) # or dos_projection='orbits'
    bs_dos_plot.get_plot(bs=bs_data, dos=dos_data)
    plt.savefig("band_dos.png")


if __name__ == "__main__":
    plot_band()
    plot_dos()
    plot_bands_dos()