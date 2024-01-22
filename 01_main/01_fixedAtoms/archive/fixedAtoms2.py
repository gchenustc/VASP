"""
#Desc: Fixing Atoms for the POSCAR format file.
#Usage: python script.py -n 16 # fixing the 16 atoms which the distance from low to high of z axis.
#Author: gchen
"""

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
    Class to fix atoms in a POSCAR file.

    Attributes:
    - input_path (str): Path to the input POSCAR file.
    - out_path (str): Path to the output POSCAR file.
    - sort_axis (str): Axis along which to sort the atoms.
    - n_fixed (int or list): Number of atoms to fix or a list of atom indices to fix.

    Methods:
    - retPoscar(): Returns the fixed POSCAR file.
    - sortAtoms(): Sorts the atoms based on the specified axis.
    - fixAtoms(): Fixes the specified number of atoms or the atoms at the specified indices.
    """

    def __init__(self, input_path='POSCAR', out_path='POSCAR_out', axis='z', n_fixed=0):
        """
        Initializes the FixAtoms object.

        Parameters:
        - input_path (str): Path to the input POSCAR file.
        - out_path (str): Path to the output POSCAR file.
        - axis (str): Axis along which to sort the atoms.
        - n_fixed (int or list): Number of atoms to fix or a list of atom indices to fix.
        """
        self.input_path = input_path
        self.out_path = out_path
        self.sort_axis = axis
        self.n_fixed = n_fixed

        self.pos_arr = np.loadtxt(input_path, dtype=np.float64, skiprows=8)
        self.pos_arr_bak = self.pos_arr.copy()
        self.const_arr = np.loadtxt(
            input_path, dtype=np.float64, skiprows=2, max_rows=3)

        # Atom information - species & count
        self.atoms_info = np.loadtxt(
            input_path, dtype='object', skiprows=5, max_rows=2).reshape((2, -1))
        self.atoms_num = self.pos_arr.shape[0]

        # Fix atoms
        self.fixAtoms()

    def retPoscar(self):
        """
        Returns the fixed POSCAR file.
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

            # Write to output file
            with open(self.out_path, 'w') as fw:
                fw.writelines(head)
                np.savetxt(fw, self.pos_arr, fmt='%s', delimiter=sep)

    def sortAtoms(self):
        """
        Sorts the atoms based on the specified axis.

        Returns:
        - sorted_indices (numpy.ndarray): Array of sorted atom indices.
        """
        # Map user input axis to numpy array axis, e.g., 'x' --> 1
        axis_range_list = ['x', 'X', 'a', 'A',
                           'y', 'Y', 'b', 'B', 'z', 'Z', 'c', 'C']
        assert self.sort_axis in axis_range_list

        axis_value_list = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
        axis_corres_dict = dict(zip(axis_range_list, axis_value_list))
        axis_value = axis_corres_dict[self.sort_axis]
        self.axis_value = axis_value

        return np.argsort(self.pos_arr_bak, 0)[:, self.axis_value]

    def fixAtoms(self):
        """
        Fixes the specified number of atoms or the atoms at the specified indices.
        """
        fixed_sign_list = []  # e.g., [['T','T','T'],['F','F','F']...]

        if isinstance(self.n_fixed, list):
            fixed_index_list = list(
                map(lambda x: x-1, self.n_fixed))  # Subtract 1 from the input numbers
        else:
            # Sort atoms and take the first n_fixed indices
            fixed_index_list = self.sortAtoms()[:self.n_fixed]

        for index, each_atom_coor in enumerate(self.pos_arr):
            if index in fixed_index_list:
                fixed_sign_list.append(['F', 'F', 'F'])
            else:
                fixed_sign_list.append(['T', 'T', 'T'])
        self.pos_arr = np.c_[self.pos_arr, np.array(
            fixed_sign_list, dtype="object")]


if __name__ == "__main__":
    for p, o, n in zip(prm.POSCAR, prm.out, prm.number):
        fix_ = FixAtoms(input_path=p, out_path=o, axis='z', n_fixed=n)
        fix_.retPoscar()
