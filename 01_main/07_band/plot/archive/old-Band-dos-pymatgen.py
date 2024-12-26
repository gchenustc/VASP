import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter,BSPlotter,BSPlotterProjected,DosPlotter
dos_vasprun=Vasprun("./vasprun.xml",)
dos_data=dos_vasprun.complete_dos
bs_vasprun=Vasprun("./vasprun.xml",parse_projected_eigen=True)
bs_data=bs_vasprun.get_band_structure(line_mode=1)


plt_1=BSDOSPlotter(bs_projection="elements",dos_projection="elements",fig_size=(20,12) )
plt_1.get_plot(bs=bs_data,dos=dos_data)
plt.savefig('band-dos.png')