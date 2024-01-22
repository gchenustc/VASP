import numpy as np
import argparse
from ase.io import write,read
"""
description:
将POSCAR转换成笛卡尔坐标
输入文件名: POSCAR
输出文件名: POSCAR_out
"""
parser = argparse.ArgumentParser(description='Conversion to Cartesian')
parser.add_argument('-c','--POSCAR', action="store", default=str('POSCAR'), dest="POSCAR",
                    help='input file name (Default=POSCAR)')
prm = parser.parse_args()

class ToCartesian(object):
    def __init__(self, input_path = 'POSCAR', out_path='POSCAR_out'):
        atoms = read(input_path)
        write('POSCAR_out', atoms, format="vasp")

if __name__=="__main__":
    ce = ToCartesian(input_path=prm.POSCAR)
