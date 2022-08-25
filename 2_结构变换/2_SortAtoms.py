import numpy as np

class SortAtoms(object):
    """
    将 POSCAR 原子按照指定轴进行排序
    SortAtoms(input_path = 'POSCAR', out_path='POSCAR_out', axis='z')
    Parameters:
    input_path: 输入文件路径
    out_path: 输出文件路径
    axis: 排序的方向，dtype=str
    
    Methods:
    retPoscar(self): 输出排序后的 POSCAR
    getTotalSortIndex(self): 得到所有原子排序的index值. eg：[3,1,0,4] --> [2,1,0,3]
    """
    def __init__(self, input_path = 'POSCAR', out_path='POSCAR_out', axis='z'):
        self.input_path = input_path
        self.out_path = out_path
        self.sort_axis = axis
        
        self.pos_arr = np.loadtxt(input_path,dtype=np.float64,skiprows=8)
        self.pos_arr_bak = self.pos_arr.copy()
        self.const_arr = np.loadtxt(input_path,dtype=np.float64,skiprows=2,max_rows=3)
        
        # 原子信息-种类&数量
        self.atoms_info = np.loadtxt(input_path,dtype='object',skiprows=5,max_rows=2).reshape((2,-1))
        self.atoms_num = self.pos_arr.shape[0]

        # 排序
        self.sortAtoms()
        
    def retPoscar(self):
        """
        输出 vasp - POSCAR 文件
        """
        with open(self.input_path,'r') as f:
            head = f.readlines()
            head_copy = head[:] # bak

            head = head[:8]
            sep = ' '*10
            head[2] = sep.join(str(self.const_arr[0].tolist()).strip('[]').split(',')) + '\n'
            head[3] = sep.join(str(self.const_arr[1].tolist()).strip('[]').split(',')) + '\n'
            head[4] = sep.join(str(self.const_arr[2].tolist()).strip('[]').split(',')) + '\n'
            head[6] = '    ' + '   '.join(str(list(map(int,self.atoms_info[1].tolist()))).strip('[]').split(',')) + '\n'  
            #head.append('Selective dynamics\n')

            # 写入
            with open(self.out_path,'w') as fw:
                fw.writelines(head)
                np.savetxt(fw,self.pos_arr,fmt='%.8f',delimiter=sep)

    def sortAtoms(self):
        # 将用户输入的排序方向与numpy数组方向对应，比如 'x'-->1
        axis_range_list = ['x','X','a','A','y','Y','b','B','z','Z','c','C']
        assert self.sort_axis in axis_range_list

        axis_value_list = [0,0,0,0,1,1,1,1,2,2,2,2]
        axis_corres_dict = dict(zip(axis_range_list,axis_value_list))
        axis_value = axis_corres_dict[self.sort_axis]
        self.axis_value = axis_value

        accu_before = 0
        accu_after = 0
        
        # 按照不同的原子类别进行排序
        for each_kind_num in map(int,self.atoms_info[1]):
            accu_after += each_kind_num
            
            pos_arr_segment = self.pos_arr[accu_before:accu_after]
            index_segment = np.argsort(pos_arr_segment,0)[:,axis_value]
            self.pos_arr[accu_before:accu_after] = pos_arr_segment[index_segment]

            accu_before += each_kind_num

        
    def getTotalSortIndex(self):
        return np.argsort(self.pos_arr_bak, 0)[:,self.axis_value]
    
                
if __name__ == "__main__":
    sort_ = SortAtoms(input_path = 'POSCAR', out_path='POSCAR_out', axis='y')
    sort_.retPoscar()
