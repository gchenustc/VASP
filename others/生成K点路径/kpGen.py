#!/usr/bin/env python
"""
A script to make kpoints file
"""

from __future__ import print_function
import seekpath
import sys
import ase
import numpy as np
import ase.io as io
import argparse
import time
# import spglib


parser = argparse.ArgumentParser(description='A script to make KPOINTS file.')
parser.add_argument("-c",
                  action="store", dest="input_file", default="POSCAR",
                  help="The crystal structure. Default: POSCAR")
parser.add_argument("-o", "--output",
                  action="store", dest="output_file", default="KPOINTS.k_path",
                  help="The file to which the new kpoints will be appended. Default: KPOINTS.kpgen.")
parser.add_argument("-r", "--resolution",
                  action="store", dest="resolution", default=0.1,
                  help="a reference target distance between neighboring k-points in the path, in units of 1/ang.")
parser.add_argument("-t", action="store_true", dest="time_reversal",
                  help="Turns off time reversal symmetry.")
parser.add_argument("-s", "--symprec", action="store", default=0.01, dest="symprec",
                  help="precision for symmetry detection [spglib]. Default: 0.01 \AA")
parser.add_argument("-e", "--explicit", action="store_false", default=True, dest="explicit",
                  help="write explicit kpoints? Default: True")
parser.add_argument("--hybrid", action="store_true", dest="hybrid",
                  help="For hybrid bandstructure calculation?")
parser.add_argument("--vdir", action="store", dest="vdir", default=None,
                  help="vacuum dir? [0->x;1->y;2->z]")
parser.add_argument("-v", action="store_true", dest="verbose", default=False,
                  help="verbose output?")
options = parser.parse_args()

# Starting
#----------------------------
if options.verbose == True:
    starttime = time.time()
    print("Starting calculation at", end='')
    print(time.strftime("%H:%M:%S on %a %d %b %Y\n"))

# read in structure
structure = io.read(options.input_file)
numbers = structure.get_atomic_numbers()
inp = (structure.cell,structure.get_scaled_positions(),numbers)

# # find symmetry with spglib
# #----------------------------
# print("Structure information [spglib]:")
# spacegroup = spglib.get_spacegroup(inp, symprec=float(options.symprec))
# print("\tSpace group number: %s" % spacegroup.split()[1])
# print("\tSpace international symbol: %s" % spacegroup.split()[0])


# Turn off time reversal symmetry if necessary
if not options.time_reversal:
    tr = True
else:
    tr = False

# get K-points
explicit_data = seekpath.get_explicit_k_path(inp,with_time_reversal=tr,
                                             reference_distance=float(options.resolution),
                                             symprec=float(options.symprec))
kpath = explicit_data['explicit_kpoints_rel']
seg = np.array(explicit_data['explicit_segments'])
labels = np.array(explicit_data['explicit_kpoints_labels'])

# return symmetry
if options.verbose == True:
    print("Structure information:")
    print("\tPrecision for finding symmetry: %8.6f \AA" % options.symprec)
    print("\tSpace group number: %s" % explicit_data['spacegroup_number'])
    print("\tSpace international symbol: %s" % explicit_data['spacegroup_international'])
    print("\nTime reversal symmtery: %s" % tr)


# 2D material?
seg_rm = []
if options.vdir!=None:
    for iseg in range(len(seg)):
        if kpath[seg[iseg,0],int(options.vdir)]!=0.0 or kpath[seg[iseg,1]-1,int(options.vdir)]!=0.0:
            seg_rm.append(iseg)
seg = np.delete(seg, seg_rm, axis=0)

# construct path
fkpath = np.array([]).reshape(0,3)
for iseg in range(len(seg)):
    fkpath=np.append(fkpath,kpath[seg[iseg,0]:seg[iseg,1]],axis=0)

# construct label path
flabels = np.array([])
for iseg in range(len(seg)):
    flabels=np.append(flabels,labels[seg[iseg,0]:seg[iseg,1]],axis=0)

if options.verbose == True:
    print("k-point path:")
    for iseg in range(len(seg)):
        print("\t%s\t(%8.6f %8.6f %8.6f)\t->\t%s\t(%8.6f %8.6f %8.6f)" %(labels[seg[iseg,0]],
                                                                         kpath[seg[iseg,0],0],
                                                                         kpath[seg[iseg,0],1],
                                                                         kpath[seg[iseg,0],2],
                                                                         labels[seg[iseg,1]-1],
                                                                         kpath[seg[iseg,1]-1,0],
                                                                         kpath[seg[iseg,1]-1,1],
                                                                         kpath[seg[iseg,1]-1,2]))

# construct label path:
#------------------------
label_path = []
for iseg in range(len(seg)):
    label_path.append(f"""{kpath[seg[iseg,0],0]:+8.6f} {kpath[seg[iseg,0],1]:+8.6f} {kpath[seg[iseg,0],2]:+8.6f} !{labels[seg[iseg,0]]}
{kpath[seg[iseg,1]-1,0]:+8.6f} {kpath[seg[iseg,1]-1,1]:+8.6f} {kpath[seg[iseg,1]-1,2]:+8.6f} !{labels[seg[iseg,1]-1]}\n\n""")

#print(''.join(label_path))

# set k-point weight to zero?
if options.hybrid:
    weight = 0.
    if options.verbose == True:
        print("For hybrid calculations.")
else:
    weight = 1./len(fkpath)

# write data
with open(options.output_file,'w') as outfile:
    if options.explicit:
        outfile.write("File generated by kp_gen.py\n")
        outfile.write(str(len(fkpath))+"\n")
        outfile.write("Reciprocal\n")
        for ipoint in range(len(fkpath)):
            outfile.write("% 8.6f % 8.6f % 8.6f %5.3f !%s\n" % (fkpath[ipoint,0],
                                                             fkpath[ipoint,1],
                                                             fkpath[ipoint,2],
                                                             weight,
                                                             flabels[ipoint]))
    else:
        outfile.write("File generated by kp_gen.py\n")
        outfile.write(str(30)+"\n")
        outfile.write("Line-mode\n")
        outfile.write("rec\n")
        outfile.write(''.join(label_path))

if options.verbose == True:
    print('Output written to {}'.format(options.output_file))
    endtime = time.time()
    runtime = endtime-starttime
    print("\nEnd of calculation.")
    print("Program was running for %.2f seconds." % runtime)

# # Write new conventions cell, to ensure compliance with the kpoints
# new_data = seekpath.get_path(inp)
#
# new_cell = ase.Atoms(positions=new_data['conv_positions'],cell=new_data['conv_lattice'],pbc=[True,True,True])
# new_cell.set_scaled_positions(new_data['conv_positions'])
# new_cell.set_atomic_numbers(new_data['conv_types'])
# numbers = new_cell.get_atomic_numbers()
#
# oup = (new_cell.cell,new_cell.get_scaled_positions(),numbers)
#
#
# # find symmetry with spglib
# #----------------------------
# print "Structure information [spglib]:"
# spacegroup = spglib.get_spacegroup(oup, symprec=float(options.symprec))
# print "\tSpace group number: %s" % spacegroup.split()[1]
# print "\tSpace international symbol: %s" % spacegroup.split()[0]
#
# if int(explicit_data['spacegroup_number'])!=int(spacegroup.split()[1][1:-1]):
#     sys.stdout.write("\033[1;31m" ) # set color red
#     print "symmetry changed!!!"
#
#
# print 'New coordinates written to CONTCAR.conventional'
# ase.io.vasp.write_vasp('CONTCAR.conventional',new_cell,sort=True,vasp5=True)
# print 'Spacegroup: {} ({})'.format(new_data['spacegroup_international'], new_data['spacegroup_number'])
# print 'Inversion symmetry?: {}'.format(new_data['has_inversion_symmetry'])
# print 'I owe you nothing.'