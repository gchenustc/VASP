# plot_vasp_band_structure.py
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter, DosPlotter

# 变量初始化
bs_vasprun = Vasprun("vasprun.xml", parse_projected_eigen=True)	
bs_data = bs_vasprun.get_band_structure(line_mode=True, force_hybrid_mode=True)  # force_hybrid_mode=True 杂化泛函开关
dos_data = bs_vasprun.complete_dos

"""
BSDOSPlotter的参数解释：cb_energy_range：导带从有带的地方开始，能量向上画到结束的区间。vb_energy_range：同理
"""
def plot_band():
    # 相当于  bands_plot = BSDOSPlotter(bs_projection='elements', dos_projection='elements')
    bands_plot = BSDOSPlotter()
    plt = bands_plot.get_plot(bs=bs_data)
    # 如果需要保存图片
    # 支持格式 png, tif, jpg....，只要改后缀名即可
    plt.savefig("bands.tif")
    plt.savefig("bands.png")

def plot_dos():
    dos_plot = DosPlotter()
    dos_plot.add_dos("Total DOS", dos_data)
    plt = dos_plot.get_plot([-4, 4], [0, 12])
    plt.savefig("dos.png")

def plot_bands_dos():
    bs_dos_plot = BSDOSPlotter(bs_projection='elements', dos_projection='elements')
    plt = bs_dos_plot.get_plot(bs=bs_data, dos=dos_data)
    plt.savefig("band_dos.png")


if __name__ == "__main__":
    plot_band()
    plot_dos()
    plot_bands_dos()