import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter, BSPlotter, BSPlotterProjected, DosPlotter

#dos_vasprun = Vasprun("./vasprun.xml")
#dos_data = dos_vasprun.complete_dos

bs_vasprun = Vasprun("./vasprun.xml", parse_projected_eigen=True)
bs_data = bs_vasprun.get_band_structure(line_mode=True)

pband_fig = BSPlotterProjected(bs=bs_data)
pband_fig = pband_fig.get_projected_plots_dots({'N':['s','p'],'Si':['s','p']})
plt.savefig('band1.png')
