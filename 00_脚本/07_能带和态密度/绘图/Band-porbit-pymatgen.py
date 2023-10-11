import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter, BSPlotter, BSPlotterProjected, DosPlotter

#dos_vasprun = Vasprun("./vasprun.xml")
#dos_data = dos_vasprun.complete_dos

bs_vasprun = Vasprun("./vasprun.xml", parse_projected_eigen=True)
bs_data = bs_vasprun.get_band_structure(line_mode=True)

plt_10=BSPlotterProjected(bs=bs_data)
plt_10.get_projected_plots_dots_patom_pmorb(selected_branches=[1,2,3,4], ylim=(-1.6,1.2),\
                                    dictio={'N':['px','py','pz'],'Si':['px','py','pz']},\
                                    dictpa={'N':[1,2,3,4,5,6,7,8,9,10,11,12],'Si':[13,14]},\
                                    sum_atoms={'N':[1,2,3,4,5,6,7,8,9,10,11,12],'Si':[13,14]})

plt.savefig('band2.png')
