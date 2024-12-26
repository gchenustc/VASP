import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun
#from matplotlib import font_manager 
#font_path = "/public/home/c1337375425/storage/fonts"
#font_files = font_manager.findSystemFonts(fontpaths=font_path)
#for file in font_files:
#    font_manager.fontManager.addfont(file)
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

vb_energy_range = 8
cb_energy_range = 8

# load data
result = Vasprun('./vasprun.xml', parse_potcar_file=False)
complete_dos = result.complete_dos
pdos_N = complete_dos.get_element_spd_dos('N')
pdos_Si = complete_dos.get_element_spd_dos('Si')

plotter = DosPlotter()
plotter.add_dos('Total DOS', complete_dos)
plotter.add_dos('N(s)', pdos_N[OrbitalType.s])
plotter.add_dos('N(p)', pdos_N[OrbitalType.p])
plotter.add_dos('Si(s)', pdos_Si[OrbitalType.s])
plotter.add_dos('Si(p)', pdos_Si[OrbitalType.p])

# ylim=(negative_value, positive_value) 可以显示spin_down 和 spin_up
plotter.get_plot(xlim=(-vb_energy_range, cb_energy_range = 8), ylim=(-12,12))
plt.savefig('dos.png')