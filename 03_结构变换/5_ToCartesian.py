import numpy as np
"""
description:
将POSCAR转换成笛卡尔坐标
输入文件名: POSCAR
输出文件名: POSCAR_out
"""

class ToCartesian(object):
    """
    转换为笛卡尔坐标，已经固定原子的POSCAR也可以使用本脚本
    """
    def __init__(self, input_path = 'POSCAR', out_path='POSCAR_out'):
        self.input_path = input_path
        self.out_path = out_path
        # 根据POSCAR第八行是否为"s"判定Selective Dynamics
        with open(input_path) as f:
            txt = f.readlines()
            if txt[7].strip().lower()[0] == 's':
                self.isfixed = True
            else:
                self.isfixed = None
                
        # 原子信息-种类&数量
        self.atoms_info = np.loadtxt(input_path,dtype='object',skiprows=5,max_rows=2).reshape((2,-1)) # np-[['N','O'],[5,5]]
        self.atoms_info[1] = self.atoms_info[1].astype('i4')
        self.atoms_num = self.atoms_info[1].sum()
        
        # 导入原子坐标
        if self.isfixed:
            self.pos_arr_ = np.loadtxt(input_path,dtype="object",skiprows=9, max_rows=self.atoms_num)
            self.pos_arr = self.pos_arr_[:,:3].astype("float64")
            self.pos_arr_suffix = self.pos_arr_[:,3:]
        else:
            self.pos_arr = np.loadtxt(input_path,dtype="float64",skiprows=8, max_rows=self.atoms_num)
#         print(self.pos_arr_suffix)

        self.pos_arr_bak = self.pos_arr.copy()
        self.const_arr = np.loadtxt(input_path,dtype=np.float64,skiprows=2,max_rows=3)
        
        # 判断原子坐标类型 - Cartesian or Direct
        if self.isfixed:
            coord_type = str(np.loadtxt('POSCAR',dtype='object',skiprows=8,max_rows=1)).lower()
        else:
            coord_type = str(np.loadtxt('POSCAR',dtype='object',skiprows=7,max_rows=1)).lower()

        # 如果是分数坐标，改成笛卡尔坐标
        if coord_type[0] == "d":
            self.pos_arr = ToCartesian.fra2car(self.pos_arr, self.const_arr)
#         print(self.pos_arr)


        
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
            if not self.isfixed:
                head[7] = 'Cartesian\n'
            else:
                self.pos_arr = np.hstack((self.pos_arr, self.pos_arr_suffix))
                head[7] = 'Selective dynamics\n'
                head.append('Cartesian\n')

            # 写入
            with open(self.out_path,'w') as fw:
                fw.writelines(head)
                np.savetxt(fw,self.pos_arr,fmt='%s',delimiter=sep)

            
if __name__=="__main__":
    ce = ToCartesian()