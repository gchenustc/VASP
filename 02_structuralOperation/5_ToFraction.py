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
        self.input_path = input_path
        self.out_path = out_path

        self.atoms_info = np.loadtxt(input_path,dtype='object',skiprows=5,max_rows=2).reshape((2,-1)) # np-[['N','O'],[5,5]]
        self.atoms_info[1] = self.atoms_info[1].astype('i4')

        atoms = read(input_path)
        self.const_arr = atoms.cell
        self.pos_arr = ToCartesian.car2fra(atoms.positions,self.const_arr)
        #write('POSCAR_out', atoms, format="vasp")
        
        # main
        self.retPoscar()
        
    @staticmethod
    def fra2car(pos_arr, const):
        """将分数坐标转换为笛卡尔坐标，传入的参数分别为分数坐标和晶格常数"""
        return np.dot(pos_arr,const)
    
    @staticmethod
    def car2fra(pos_arr, const):
        """将笛卡尔坐标转换为分数坐标，传入的参数分别为笛卡尔坐标和晶格常数"""
        return np.dot(pos_arr, np.linalg.inv(const))

    def retPoscar(self):
        """
        输出 vasp - POSCAR 文件
        """
        with open(self.input_path,'r') as f:
            head = f.readlines()
            head_copy = head[:] # bak

            head = head[:8]
            sep = ' '*10
            head[2] = sep.join(map(str,self.const_arr[0].tolist())) + '\n'
            head[3] = sep.join(map(str,self.const_arr[1].tolist())) + '\n'
            head[4] = sep.join(map(str,self.const_arr[2].tolist())) + '\n'
            head[5] = '    ' + ' '.join(self.atoms_info[0].tolist()) + '\n'
            head[6] = '    ' + ' '.join(list(map(str,self.atoms_info[1].tolist()))) + '\n'
            head[7] = 'Direct\n'

            # 写入
            with open(self.out_path,'w') as fw:
                fw.writelines(head)
                np.savetxt(fw,self.pos_arr,fmt='%s',delimiter=sep)

            
if __name__=="__main__":
    ce = ToCartesian(input_path=prm.POSCAR)
