import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter, BSPlotter, BSPlotterProjected, DosPlotter
#from matplotlib import font_manager 
#font_path = "/public/home/c1337375425/storage/fonts"
#font_files = font_manager.findSystemFonts(fontpaths=font_path)
#for file in font_files:
#    font_manager.fontManager.addfont(file)
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

bs_vasprun = Vasprun("./vasprun.xml", parse_projected_eigen=True)
bs_data = bs_vasprun.get_band_structure(line_mode=True) # 杂化泛函需要同时打开：line_mode=True, force_hybrid_mode=True，需要注意的是杂化泛函的能带计算也要读入k点信息，所以保证在KPOINTS文件里的高对称k点的后面添加标记，这样get_band_structure函数可以自动的识别这些高对称点，否则，需要手动载入包含k点信息的字典。

plt_10=BSPlotterProjected(bs=bs_data)
plt_10.get_projected_plots_dots_patom_pmorb(selected_branches=[1,2,3,4], ylim=(-1.6,1.2),\
                                    dictio={'N':['px','py','pz'],'Si':['px','py','pz']},\
                                    dictpa={'N':[1,2,3,4,5,6,7,8,9,10,11,12],'Si':[13,14]},\
                                    sum_atoms={'N':[1,2,3,4,5,6,7,8,9,10,11,12],'Si':[13,14]})

plt.savefig('band2.png')
