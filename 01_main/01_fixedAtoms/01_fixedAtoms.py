DESC=USAGE=\
"""
Author: gchen
Time: 20240120
Version: 3.0
Description:
    This script allows you to fix specific atoms in a crystal structure by applying constraints.
    You can customize the fixed atoms based on different criteria such as element, number, or ratio.
Example usage:
    python script.py -p 24                                     # fixing 24 atoms from low to high position along z axis
    python script.py -p 24 -a y                                # fixing 24 atoms from low to high position along y axis
    python script.py -p 24 -o POSCAR                           # fixing 24 atoms and assign the output file name
    python script.py -s xxx.vasp -p 24 -o POSCAR               # fixing 24 atoms and assign the input file name and output file name
    python script.py -p 24 -r                                  # fixing 24 atoms from top of z-axis
    python script.py -p 24 -a y -r                             # fixing 24 atoms from top of y-axis
    python script.py -p 24 -g                                  # fixing 24 atoms, but reserve the original fixed way.
    python script.py -s POSCAR -d 1 0 1 -p 24                  # fixing 24 atoms and assign direction of "1 0 1"
    python script.py -s POSCAR -d 1 0 1 -t 2 -p "N"            # fixing all N atoms and assign direction of "1 0 1"
    python script.py -s POSCAR -d 1 0 1 -t 3 -p 0.5 -q 20      # fixing 50% atoms out of 20 atoms. 
    python script.py -s POSCAR -o xxx.vasp -a y -d 1 0 1 -t 3 -p 0.5 -q 20 -r -g      
                                                               # fixing 50% atoms out of 20 atoms from top of y axis for the fixed way of "F T F" and reserve the original fixed way, the input file is POSCAR, the output file is xxx.vasp. 
Arguments:
    -s, --stru_path  
                            Path to the structure file (e.g., POSCAR). Default: POSCAR
    -o, --output_path
                            The output file name. Default: POSCAR
    -d, --fixed_direction
                            Fixed direction as a list of integers. Default: 1 1 1
    -a  --axis                           
                            The fixed axis. Default: z. choices: {x(X), y(Y), z(Z)}
    -t, --fixed_type      
                            Type of fixing:
                            1: Fixing according to number
                            2: Fixing according to element
                            3: Fixing according to ratio. Default: 1
                            choices: {1, 2, 3}
    -p, --param1        
                            Parameter 1:
                            If fixed_type = 1, param1 is the fixed number
                            If fixed_type = 2, param1 is the element name
                            If fixed_type = 3, param1 is the fixed ratio. Default: "0"
    -q, --param2         
                            Parameter 2: 
                            Only used when fixed_type = 3. It is the total number of atoms to calculate the fixed ratio. Default: 0
    -r
                            If reversed z axis (from top to bottom of z), only used when fixed_type=1/3. Default: None
    -g
                            If reserve the fixed atoms for original file (Noted, if original file is not fixed, this script is not fixed it, it will fix it. But the if original file is fixed, the new is not fixed, it will also fix following the original fixed way.). Default: None
"""

import argparse
import numpy as np
from ase.io import read
from ase.io import write
import math
from ase.constraints import FixAtoms, FixedLine, FixedPlane


def judgeFixedDirection(direction):
    n = sum(i for i in direction if i != 0)
    if n == 1:
        return FixedPlane, list(map(lambda x: 0 if x == 0 else 1, direction))
    elif n == 2:
        return FixedLine, list(map(lambda x: 0 if x != 0 else 1, direction))
    else:
        return FixAtoms, [1, 1, 1]


def fixedOperation(stru, indices, direction=None, reserve_old=False):
    if direction is None:
        direction = [1, 1, 1]
    constraint = []
    
    if reserve_old:
        constraint.extend(stru.constraints)

    func, dire = judgeFixedDirection(direction)
    for idx in indices:
        # if the direction=[1,1,1], the func == FixAtoms and don't need have direction parameter.
        tmp = func([idx]) if dire == [1, 1, 1] else func(idx, dire)
        constraint.append(tmp)
        
    stru.set_constraint(constraint)


def fixAccordingElement(stru, direction=None, ele=None, reserve_old=False):
    if direction is None:
        direction = [1, 1, 1]
    indices = [atom.index for atom in stru if atom.symbol == ele]
    fixedOperation(stru, indices, direction=direction, reserve_old=reserve_old)


def fixAccordingNumber(stru, direction=None, n=0, reverse=0, reserve_old=False, axis="z"):
    if direction is None:
        direction = [1, 1, 1]
    
    assert axis in ["x", "y", "z", "X", "Y", "Z"]
    if axis.lower() == "z":
        axis_value = 2
    elif axis.lower() == "x":
        axis_value = 0
    elif axis.lower() == "y":
        axis_value = 1

    if reverse:
        ranked_index = stru.get_positions()[:, axis_value].argsort()[::-1]
    else:
        ranked_index = stru.get_positions()[:, axis_value].argsort()
    indices = ranked_index[:n]

    fixedOperation(stru, indices, direction=direction, reserve_old=reserve_old)


def fixAccordingRatio(stru, direction=None, ratio=0, n_total=0, reverse=0, reserve_old=False, axis="z"):
    if direction is None:
        direction = [1, 1, 1]
        
    assert axis in ["x", "y", "z", "X", "Y", "Z"]
    if axis.lower() == "z":
        axis_value = 2
    elif axis.lower() == "x":
        axis_value = 0
    elif axis.lower() == "y":
        axis_value = 1
        
    if reverse:
        ranked_index = stru.get_positions()[:, axis_value].argsort()[::-1]
    else:
        ranked_index = stru.get_positions()[:, axis_value].argsort()

    # the farthest atoms' z value out of total.
    indices = ranked_index[:math.floor(n_total * ratio)]
    fixedOperation(stru, indices, direction=direction, reserve_old=reserve_old)


if __name__ == "__main__":
    # Create ArgumentParser object and add arguments
    parser = argparse.ArgumentParser(description=DESC, usage=USAGE, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-s", "--stru_path", help="Path to the structure file (e.g., POSCAR). Default: POSCAR", default="POSCAR")
    parser.add_argument("-d", "--fixed_direction", nargs="+", type=int, help="Fixed direction as a list of integers. Default: 1 1 1", default=[1, 1, 1])
    parser.add_argument("-a", "--axis", type=str, choices=["x", "X", "y", "Y", "z", "Z"], help="The fixed axis. Default: z", default="z")
    parser.add_argument("-t", "--fixed_type", type=int, choices=[1, 2, 3], help="Type of fixing:\n1: Fixing according to element\n2: Fixing according to number\n3: Fixing according to ratio. Default: 1", default=1)
    parser.add_argument("-p", "--param1", type=str, help="Parameter 1:\nIf fixed_type = 1, param1 is the element name\nIf fixed_type = 2, param1 is the fixed number\nIf fixedtype = 3, param1 is the fixed ratio. Default: 0", default="0")
    parser.add_argument("-q", "--param2", type=int, help="Parameter 2: Only used when fixed_type = 3. It is the total number of atoms to calculate the fixed ratio. Default: 0", default=0)
    parser.add_argument("-r", action='store_const', const=True, help="If reversed z axis (from top to bottom of z), only used when fixed_type=1/3. Default: None")
    parser.add_argument("-g", action='store_const', const=True, help="If reserve the fixed atoms for original file (Noted, if original file is not fixed, this script is not fixed it, it will fix it. But the if original file is fixed, the new is not fixed, it will also fix following the original fixed way.). Default: None")
    parser.add_argument("-o", "--output_path", type=str, help="The output file name. Default: POSCAR", default="POSCAR")
   
    args = parser.parse_args()

    stru = read(args.stru_path)
    if args.fixed_type == 1:
        fixAccordingNumber(stru, direction=args.fixed_direction, n=int(args.param1), reverse=args.r, reserve_old=args.g, axis=args.axis)
    elif args.fixed_type == 2:
        fixAccordingElement(stru, direction=args.fixed_direction, ele=args.param1, reserve_old=args.g)
    elif args.fixed_type == 3:
        fixAccordingRatio(stru, direction=args.fixed_direction, ratio=float(args.param1), n_total=args.param2, reverse=args.r, reserve_old=args.g, axis=args.axis)
    write(args.output_path, stru, format="vasp")