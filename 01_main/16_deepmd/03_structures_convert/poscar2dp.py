from dpdata import LabeledSystem,MultiSystems
from glob import glob
"""
DESC: convert the POSCARs/OUTCARs set into deepmd's data. It can Process multi systems. Please note to modify the vasp_path and dp_path in this script.
USAGE:
        python /path/to/script.py
"""
vasp_path = './*/*/*/OUTCAR'
dp_path = "./deepmd"

fs=glob(vasp_path)  # remeber to change here !!!
ms=MultiSystems()
for f in fs:
    try:
        ls=LabeledSystem(f)
    except:
        print(f)
    if len(ls)>0:
        ms.append(ls)

ms.to_deepmd_raw(dp_path)
ms.to_deepmd_npy(dp_path)