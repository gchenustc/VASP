import numpy as np
from itertools import permutations # A_{n,n}
from itertools import combinations # C_{n,n}
from itertools import combinations_with_replacement # 可放回

class HollowCuboid(object):
    """
    挖出 POSCAR 指定位置的 长方体
    使用方法:
    1. 指定输入输出文件路径
    hc = HollowCuboid(input_path = 'POSCAR', out_path='POSCAR_out')
    2. 指定需要挖出的长方体中心以及长方体的晶格常数，也可使长方体旋转,rotate_axis类型使角度，rotate_angle是rad角，可以用np.pi/n
    hc.operate(center_point=[5,5,5], lattice_constant=[10,10,10], rotate_axis=[0,0,1], rotate_angle=np.pi/4)
    3. 输出
    hc.retPoscar()
    """
    def __init__(self, input_path = 'POSCAR', out_path='POSCAR_out'):
        self.input_path = input_path
        self.out_path = out_path
        
        self.pos_arr = np.loadtxt(input_path,dtype=np.float64,skiprows=8)
        self.pos_arr_bak = self.pos_arr.copy()
        self.const_arr = np.loadtxt(input_path,dtype=np.float64,skiprows=2,max_rows=3)
        
        # 判断原子坐标类型 - Cartesian or Direct
        coord_type = str(np.loadtxt('POSCAR',dtype='object',skiprows=7,max_rows=1)).lower()
        # 如果是分数坐标，改成笛卡尔坐标
        if coord_type[0] == "d":
            self.pos_arr = HollowCuboid.fra2car(self.pos_arr, self.const_arr)
        
        # 原子信息-种类&数量
        self.atoms_info = np.loadtxt(input_path,dtype='object',skiprows=5,max_rows=2).reshape((2,-1)) # np-[['N','O'],[5,5]]
        self.atoms_num = self.pos_arr.shape[0]

    @staticmethod
    def fra2car(pos_arr, const):
        """将分数坐标转换为笛卡尔坐标，传入的参数分别为分数坐标和晶格常数"""
        return np.dot(pos_arr,const)
    
    @staticmethod
    def car2fra(pos_arr, const):
        """将笛卡尔坐标转换为分数坐标，传入的参数分别为笛卡尔坐标和晶格常数"""
        return np.dot(pos_arr, np.linalg.inv(const))

    @staticmethod
    def periodBoundary(pos_arr, const_arr, direction=None):
        """
        返回扩一倍胞后的等价位置
        direction 是扩胞后周期性边界条件的方向 eg: [[1,1,0],[1,0,1]...]
        假如 pos_arr.shape = (10,3)
        则返回数组的shape: ret.shape = (10,27,3), 27代表每个原子等价的27个位置
        const_arr: 晶格常数数组, const_arr.shape=(3,3)
        """
        ret = pos_arr[:,np.newaxis,:] # 初始 ret.shape = [10,1,3]，方便之后 concatenate

        direct = [] # [(1,1,0),(1-1,0)...]，(1,1,0)代表沿着x，y轴正方向扩胞的
        for i in combinations_with_replacement([1,-1,0],3): # 有放回地从[1,-1,0]中取三个值
            for j in permutations(i,3): # A_{3,3}排列组合
                direct.append(j)
        direct = list(set(direct)) # remove repeat elements
        direct.remove((0,0,0)) # remove self

        for dir_ in direct: # (1,1,0)
            dir_arr = np.array(dir_)

            const_sign_change = const_arr*(dir_arr.reshape((3,-1)))# 改变const_arr 的符号
            dir_arr_abs = np.apply_along_axis(lambda x: abs(x), arr=dir_arr,axis=0) # 将 取dir_内元素的绝对值
            const_index = np.argwhere(dir_arr_abs==1).flatten() 
            vector_arr = const_sign_change[const_index] # shape=n*3，原始胞内的每个原子坐标依次加上vector_arr中的坐标即是其他周期的等价坐标

            new_cell = pos_arr.copy()
            for vector in vector_arr:
                new_cell += vector
            new_cell = new_cell[:,np.newaxis,:] # 方便concatenate 的操作

            ret = np.concatenate((ret,new_cell),1)
        return ret


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
    def isInBox(coord,x,y,z):
        """coord 是坐标 np.array([a,b,c]),当a在x[0]与x[1]之间，同时b在y[0]和y[1]之间，c在...., 时返回True"""
        lst = [x,y,z]
        for index,i in enumerate(coord):
            if not lst[index][0] <= i <= lst[index][1]:
                return False
        return True

    @staticmethod
    def getNotInCubicIndex(pos_arr, center_point, lattice_constant):
        """
        返回晶胞中不在以center_point为中心的长方体中的原子索引: eg: [0,0,1,1,0,0],value = 1 的原子是不在长方体内的原子索引
        pos_arr 可是是二维 shape=(n,3) 的形式，也可是三维 shape=(n,m,3)，其中m是等价的位置，比如(1,m,3)，表示索引为1的原子有m个等价的情况
    
        """
        if pos_arr.ndim == 2:
            pos_arr = pos_arr[:,np.newaxis,:]
    
        ret = np.ones(pos_arr.shape[0]) # one-hot 形式
        x_range = center_point[0] - lattice_constant[0]/float(2), center_point[0] + lattice_constant[0]/float(2)
        y_range = center_point[1] - lattice_constant[1]/float(2), center_point[1] + lattice_constant[1]/float(2)
        z_range = center_point[2] - lattice_constant[2]/float(2), center_point[2] + lattice_constant[2]/float(2)

        for index,i in enumerate(pos_arr):
            for each_coord in i:
                if HollowCuboid.isInBox(each_coord,x_range,y_range,z_range):
                    ret[index] = 0
        return ret
    
    @staticmethod
    def normalize(arr):
        "归一化一维矢量,可以传入list"
        norm = np.linalg.norm(np.array(arr))
        return np.array(arr/norm)
    
    @staticmethod
    def getRotationMatrix(b, theta):
        """
        由旋转轴 b 和旋转角度 theta 求旋转矩阵 R
        theta 单位: rad
        """
        b_vect = HollowCuboid.normalize(b)
        c = np.cos(theta)
        s = np.sin(theta)
        a = 1 - np.cos(theta)
        x,y,z = b_vect

        Rx = np.array([a * x ** 2 + c, a * y * x -s * z, a * z * x + s * y])
        Ry = np.array([a * x * y + s * z, a * y ** 2 + c, a * z * y - s * x])
        Rz = np.array([a * x * z - s * y, a * y * z + s * x, a * z ** 2 + c])

        return np.array([Rx, Ry, Rz])

    def operate(self, center_point, lattice_constant, rotate_axis=None, rotate_angle=0):
        """
        lattice_constant: 小立方体的晶格常数
        """
        center_point = np.array(center_point)
        pos_arr_expand = HollowCuboid.periodBoundary(self.pos_arr, self.const_arr, direction=None).reshape((-1,3))
        if rotate_axis:
            # 得到大胞的旋转矩阵(因为是相对于小胞旋转，故旋转角度取负)
            rotate_matrix = HollowCuboid.getRotationMatrix(rotate_axis, -rotate_angle)
            # 得到旋转后的大胞
            pos_arr_rotate = np.dot(pos_arr_expand,rotate_matrix)
            # 调整大胞位置-从自身沿着某轴旋转到沿着小胞中心的某个轴旋转
                # 旋转后的小胞中心点
            center_point_ = np.dot(center_point, rotate_matrix)
            pos_arr_rotate += center_point - center_point_
            notinbox_ont_hot = HollowCuboid.getNotInCubicIndex(pos_arr_rotate.reshape((-1,27,3)), center_point, lattice_constant)
        else:
            notinbox_ont_hot = HollowCuboid.getNotInCubicIndex(pos_arr_expand.reshape((-1,27,3)), center_point, lattice_constant)
        
        # 改变 POSCAR 信息
        notinbox_index = np.argwhere(notinbox_ont_hot==1).flatten()
        self.pos_arr = self.pos_arr[notinbox_index] # 改变坐标
        
            # 改变 POSCAR 第六行的原子数量信息
        n_atoms_cumsum = np.cumsum(self.atoms_info[1].astype('i4'))
        n_atoms_per_kind = self.atoms_info[1].astype('i4')
        for index,i in enumerate(notinbox_ont_hot):
            if i==0:
                for index_,n_each_kind in enumerate(n_atoms_cumsum):
                    if index_ == 0:
                        if index < n_each_kind:
                            n_atoms_per_kind[0] -= 1
                    else:
                        if n_atoms_cumsum[index_-1] <= index < n_each_kind:
                            n_atoms_per_kind[index_] -= 1
        self.atoms_info[1] = n_atoms_per_kind
    
if __name__=="__main__":
    hc = HollowCuboid()
    hc.operate(center_point=[0,0,0], lattice_constant=[10,10,10], rotate_axis=[0,0,1], rotate_angle=np.pi/4)
    hc.retPoscar()
