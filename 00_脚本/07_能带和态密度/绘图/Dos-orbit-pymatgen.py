import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun

# load data
result = Vasprun('./vasprun.xml', parse_potcar_file=False)
complete_dos = result.complete_dos
pdos_N = complete_dos.get_element_spd_dos('N')
pdos_Si = complete_dos.get_element_spd_dos('Si')

plotter = DosPlotter()
plotter.add_dos('Total DOS', result.tdos)
plotter.add_dos('N(s)', pdos_N[OrbitalType.s])
plotter.add_dos('N(p)', pdos_N[OrbitalType.p])
plotter.add_dos('Si(s)', pdos_Si[OrbitalType.s])
plotter.add_dos('Si(p)', pdos_Si[OrbitalType.p])


plotter.get_plot(xlim=(-10,5), ylim=(-10,10))
plt.savefig('dos.png')