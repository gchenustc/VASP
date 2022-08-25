import numpy as np

class CellExpansion(object):
    """
    扩胞
    CellExpansion(input_path = 'POSCAR', out_path='POSCAR_out', cell = [2,2,2])
    """
    def __init__(self, input_path = 'POSCAR', out_path='POSCAR_out', cell = [1,1,1]):
        self.input_path = input_path
        self.out_path = out_path
        self.cell =cell
        
        self.pos_arr = np.loadtxt(input_path,dtype=np.float64,skiprows=8)
        self.pos_arr_bak = self.pos_arr.copy()
        self.const_arr = np.loadtxt(input_path,dtype=np.float64,skiprows=2,max_rows=3)
        
        # 判断原子坐标类型 - Cartesian or Direct
        coord_type = str(np.loadtxt('POSCAR',dtype='object',skiprows=7,max_rows=1)).lower()
        # 如果是分数坐标，改成笛卡尔坐标
        if coord_type[0] == "d":
            self.pos_arr = CellExpansion.fra2car(self.pos_arr, self.const_arr)
        
        # 原子信息-种类&数量
        self.atoms_info = np.loadtxt(input_path,dtype='object',skiprows=5,max_rows=2).reshape((2,-1)) # np-[['N','O'],[5,5]]
        self.atoms_info[1] = self.atoms_info[1].astype('i4')
        self.atoms_num = self.pos_arr.shape[0]
        
        # main
        self.expandCell()
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
            head[7] = 'Cartesian\n'
            #head.append('Selective dynamics\n')

            # 写入
            with open(self.out_path,'w') as fw:
                fw.writelines(head)
                np.savetxt(fw,self.pos_arr,fmt='%.8f',delimiter=sep)
                
    @staticmethod
    def expandAlongSingleAxis(pos_arr, const_arr, direction=[1,0,0], value=3):
        """
        """
        pos = pos_arr
        ret = pos_arr[:,np.newaxis,:].copy() # 初始 ret.shape = [10,1,3]，方便之后 concatenate
        direct_index = np.argwhere(np.array(direction)==1)[0][0] # 返回0,1或者2
        vector = const_arr[direct_index]
        
        for i in range(value):
            pos_arr = pos_arr + vector  # numpy 中 a = a+b 和 a +=b 不一样，前者a的地址更新，后者不会
            pos_arr_concat = pos_arr[:,np.newaxis,:]
            ret = np.concatenate((ret,pos_arr_concat),1)
            #print(pos[:10], pos_arr[:10], ret[:10])
        return ret
    
    def expandCell(self):
        for index,i in enumerate(self.cell):
            direction = np.zeros(3)
            direction[index] = 1
            self.pos_arr = CellExpansion.expandAlongSingleAxis(self.pos_arr, self.const_arr, direction, i-1).reshape((-1,3))
        
        # 原子数目改变
        self.atoms_info[1] *=  np.multiply.accumulate(np.array(self.cell))[2]
        # 胞长改变
        self.const_arr *= np.array(self.cell).reshape((3,1))
            
if __name__=="__main__":
    ce = CellExpansion(cell=[2,2,2])