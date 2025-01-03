from ase.io import read, write
import ase.db
import sys
"""
author: gchen
time: 20241212
desc: export structure from fps.db by id, and save it as POSCAR file, only one argument (one id) is needed
usage: python export_fp_stru.py id
"""

if __name__ == "__main__":
    if len(sys.argv) > 2:
        print("only one argument is needed!")
        sys.exit(1)
    if not sys.argv[1].isdigit():
        print("argument must be a number!")
        sys.exit(1)
    
    id = int(sys.argv[1])
    mydb = ase.db.connect("fps.db")
    atoms = mydb.get(id).toatoms()
    write("POSCAR", atoms, format="vasp")



