import numpy as np
import argparse

# argparser
parser = argparse.ArgumentParser(description="None")
parser.add_argument(
    '-c', '--POSCAR', default=["POSCAR"], nargs="+", help="Specify the POSCAR file")
parser.add_argument(
    '-o', '--out', default=["POSCAR"], nargs="+", help="Specify the out file")
parser.add_argument(
    '-n', '--number', default=[12], type=int, nargs="+", help="Specify fixed number of atoms")
prm = parser.parse_args()


class FixAtoms(object):
    """
    将 POSCAR 原子进行固定
    SortAtoms(input_path = 'POSCAR', out_path='POSCAR_out', axis='z',  n_fixed = 0)
    Parameters:
    input_path: 输入文件路径
    out_path: 输出文件路径
    axis: 排序的方向，dtype=str
    n_fixed: 按照排序的方向固定原子的数目，n_fixed的类型如果是list，则不进行排序，固定列表中的对应序号的原子

    Methods:
    retPoscar(self): 输出固定后的 POSCAR
    """

    def __init__(self, input_path='POSCAR', out_path='POSCAR_out', axis='z', n_fixed=0):
        self.input_path = input_path
        self.out_path = out_path
        self.sort_axis = axis
        self.n_fixed = n_fixed

        self.pos_arr = np.loadtxt(input_path, dtype=np.float64, skiprows=8)
        self.pos_arr_bak = self.pos_arr.copy()
        self.const_arr = np.loadtxt(
            input_path, dtype=np.float64, skiprows=2, max_rows=3)

        # 原子信息-种类&数量
        self.atoms_info = np.loadtxt(
            input_path, dtype='object', skiprows=5, max_rows=2).reshape((2, -1))
        self.atoms_num = self.pos_arr.shape[0]

        # 固定
        self.fixAtoms()

    def retPoscar(self):
        """
        输出 vasp - POSCAR 文件
        """
        with open(self.input_path, 'r') as f:
            head = f.readlines()
            head_copy = head[:]  # bak

            head = head[:8]
            sep = ' '*2
            head[2] = sep.join(map(str, self.const_arr[0].tolist())) + '\n'
            head[3] = sep.join(map(str, self.const_arr[1].tolist())) + '\n'
            head[4] = sep.join(map(str, self.const_arr[2].tolist())) + '\n'
            head[6] = '    ' + ' '.join(self.atoms_info[1].tolist()) + '\n'
            head.insert(7, 'Selective dynamics\n')

            # 写入
            with open(self.out_path, 'w') as fw:
                fw.writelines(head)
                np.savetxt(fw, self.pos_arr, fmt='%s', delimiter=sep)

    def sortAtoms(self):
        # 将用户输入的排序方向与numpy数组方向对应，比如 'x'-->1
        axis_range_list = ['x', 'X', 'a', 'A',
                           'y', 'Y', 'b', 'B', 'z', 'Z', 'c', 'C']
        assert self.sort_axis in axis_range_list

        axis_value_list = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
        axis_corres_dict = dict(zip(axis_range_list, axis_value_list))
        axis_value = axis_corres_dict[self.sort_axis]
        self.axis_value = axis_value

        return np.argsort(self.pos_arr_bak, 0)[:, self.axis_value]

    def fixAtoms(self):
        fixed_sign_list = []  # 比如[['T','T','T'],['F','F','F']...]

        if isinstance(self.n_fixed, list):
            fixed_index_list = list(
                map(lambda x: x-1, self.n_fixed))  # 输入的数字要减一
        else:
            # 对原子进行排序后取前 n_fixed 个索引
            fixed_index_list = self.sortAtoms()[:self.n_fixed]

        for index, each_atom_coor in enumerate(self.pos_arr):
            if index in fixed_index_list:
                fixed_sign_list.append(['F', 'F', 'T'])
            else:
                fixed_sign_list.append(['T', 'T', 'T'])
        self.pos_arr = np.c_[self.pos_arr, np.array(
            fixed_sign_list, dtype="object")]


if __name__ == "__main__":
    for p, o, n in zip(prm.POSCAR, prm.out, prm.number):
        fix_ = FixAtoms(input_path = p, out_path=o, axis='z', n_fixed=n)
        fix_.retPoscar()
