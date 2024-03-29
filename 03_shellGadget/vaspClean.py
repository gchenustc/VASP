#!/usr/bin/env python
"""
Desc: A script to clean current working dir.
Useage:
vasp_clean.py
vasp_clean.py -f   # force clean and no remind
vasp_clean.py -a REPORT XDATCAR  # reserve REPORT and XDATCAR
vasp_clean.py -f -a REPORT XDATCAR 
"""

from __future__ import print_function
import os
import sys
import argparse

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to clean current working dir.')
parser.add_argument('-f','--force', action="store_true", default=False, dest='force',
                    help='Remove without warning? (Default=False)')
parser.add_argument('-a','--add', nargs='+', action="store", default=None, dest="additional",
                    help='Additional files to retain.')

prm = parser.parse_args()

# get a list of objects under current dir
dir_file = os.listdir('./')

# target files
files = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'vdw_kernel.bindat', '.pbs', 'wanier90', '.sh', '.py']

if prm.additional != None:
    files += prm.additional

# exclude target files
for i in files:
    dir_file = [x for x in dir_file if not x.startswith(i)]
    dir_file = [x for x in dir_file if not x.endswith(i)]

# exclude directories
dir_file = [i for i in dir_file if not os.path.isdir(i)]

if prm.force == False:
    # Warning is very much needed!
    print('Removing:\n\t','\n\t'.join(dir_file))
    if float(sys.version.split()[0][:3]) < 3.0:
        val = raw_input("Are you sure? (y/n) ")
    else:
        val = input("Are you sure? (y/n) ")
    if val == 'y':
        # remove files
        for file in dir_file:
            os.remove(file)
        print('done.')
    elif val == 'n':
        print("exiting...")
    else:
        print('Input not accepted, exiting...')
        exit()
else:
    # remove files
    for file in dir_file:
        os.remove(file)
