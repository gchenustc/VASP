"""
Description:
    The provided code is a program that generates doped structures based on a template structure. It replaces specified elements in the template structure with a doping element to create new structures. The program ensures that the generated structures are not duplicated with previously generated structures or existing structures in a given directory.

Usage:
    To use the code, follow these steps:
    Make sure you have Python installed on your system.
    Save the all code (include generalFun.py package dependent) in Python files (e.g., generalFun.py and doped_structure_generator.py).
    Open a command prompt or terminal and navigate to the directory where the Python file is saved.
    Run the following command to execute the code:
    
        python doped_structure_generator.py -i <input_file> -de <doped_element> dp <doping_element> -n <n_doped> -o <n_out> [-s <index_start>] [-p <prefix>] [-l <existed_strus_dir_name>]
        
    Replace the arguments within <angle brackets> with the appropriate values:

        <input_file>: The path to the input structure file. default: POSCAR
        <doped_element>: The element to be doped. 
        <doping_element>: Element that replacing doped_element.
        <n_doped>: The number of doping elements to replace in the structure. default: 1
        <n_out>: The desired number of output doped structures. default: 1
        <index_start> (optional): The starting index for the output structure filenames. Default: 1.
        <prefix> (optional): The prefix to be used for the output structure filenames. Default: "doped_structure_".
        <existed_strus_dir_name> (optional): The directory name containing existing structures. Default: "existed_structures".
    
    Note: 
        The input structure file should be in the VASP format.

        The program will generate the specified number of doped structures based on the provided inputs. The generated structures will be saved as POSCAR files in the current directory.

        Make sure to replace the arguments with appropriate values according to your specific use case.

"""

from pymatgen.core import Structure as Sr
import generalFun as gf 
from itertools import combinations
import os
import random
import time
import logging
import argparse

# Configure logging to write messages to the log.txt file
logging.basicConfig(level=logging.INFO, filename="log.txt", filemode="a",
                    format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')

def get_doped_elements_index(tem_stru, doped_element):
    """Given the element to be replaced, return the indices of the element in the structure."""
    return [
        index
        for index, stru in enumerate(tem_stru)
        if stru.species_string == doped_element
    ]


def doping(tem_stru, doping_element, doped_index):
    """Generate a doped structure based on the given structure, doping element, and indices to be doped."""
    stru_doped = tem_stru.copy()  # Create a copy of the structure
    for index in doped_index:
        stru_doped[index] = doping_element  # Replace the specified indices with the doping element
    stru_doped.sort()  # Sort the structure
    return stru_doped


def duplicated_check(stru, strus_list):
    """Check for duplicated structures. Return False if no duplicated structures are found."""
    if len(strus_list) == 0:
        return False

    return any(gf.symmetryCheck(stru_, stru) for stru_ in strus_list)


def get_old_sturs_path(existed_strus_dir_name):
    """Get the paths of the structures in the subdirectories of existed_strus_dir_name. The generated structures will not be equivalent to the structures in the subdirectories."""
    if not os.path.exists(existed_strus_dir_name):
        return []  # Return an empty list if the directory does not exist
    old_strus_dirs_inner_list = list(os.walk(existed_strus_dir_name))[0][1]  # Get the subdirectories inside the directory
    existed_strus_path_dirs = [os.path.join(existed_strus_dir_name, i) for i in old_strus_dirs_inner_list]  # Get the paths of the subdirectories
    existed_strus_path_list_files = [list(os.walk(i))[0][2] for i in existed_strus_path_dirs]  # Get the files in the subdirectories
    old_strus_list = []
    for dir_, files in zip(existed_strus_path_dirs, existed_strus_path_list_files):
        old_strus_list.extend(os.path.join(dir_, file) for file in files)
    return old_strus_list


def get_all_doping_strus(tem_stru, doped_element, doping_element, n_doped, n_out, index_start, prefix, existed_strus_dir_name):
    """
    Generate all doped structures, excluding duplicated structures and structures in existed_strus_dir_name.
    tem_stru: Template structure for doping.
    doped_element: Element to be doped.
    doping_element: Element that replacing doped_element.
    n_doped: Number of doping elements.
    n_out: Number of output structures.
    index_start: Starting index of the output structure filenames.
    prefix: Prefix of the output structure filenames.
    existed_strus_dir_name: Structures in the subdirectories of existed_strus_dir_name will not be equivalent to the generated structures.
    """
    assert n_out > 0  # Ensure that the number of output structures is positive
    doped_elements_index = get_doped_elements_index(tem_stru, doped_element)

    assert len(doped_elements_index) >= n_doped  # Ensure that the number of elements to be replaced is greater than or equal to the number of doping elements

    old_strus_list = get_old_sturs_path(existed_strus_dir_name)
	
    existed_strus_list = [Sr.from_file(stru_path) for stru_path in old_strus_list]  # Load the existing structures from the file paths

    all_doping_strus = []  # Initialize a list to store all the doped structures
    count = 0  # Initialize a counter for the number of generated structures

    while count < n_out:
        random.shuffle(doped_elements_index)
        doped_index = list(combinations(doped_elements_index, n_doped))[0]  # Generate one combinations of indices to be doped
        doped_stru = doping(tem_stru, doping_element, doped_index)  # Generate a doped structure based on the indices   
        
        logging.info(f"checking no.{count+1} stru.")
        if not duplicated_check(doped_stru, all_doping_strus) and not duplicated_check(doped_stru, existed_strus_list):
            logging.info("\t... This is the required structure ...")
            # Check if the doped structure is not duplicated with previously generated structures or existing structures
            stru_name = prefix + str(index_start + count) + ".vasp"  # Generate the filename for the output structure
            doped_stru.to(filename=stru_name,fmt="POSCAR")
            all_doping_strus.append(doped_stru)  # Add the doped structure to the list
            
            count += 1  # Increment the counter
            if count >= n_out:
                break  # Break the loop if the desired number of output structures has been reached
        else:
            logging.info("\t... repeated ...")

    return all_doping_strus  # Return the list of all doped structures

def main():
    parser = argparse.ArgumentParser(description='Generate doped structures')
    parser.add_argument('-i', '--input', type=str, default="POSCAR", help='Input structure file')
    parser.add_argument('-de', '--doped_element', type=str, required=True, help='Element to be doped')
    parser.add_argument('-dp', '--doping_element', type=str, required=True, help='Element to be doped')
    parser.add_argument('-n', '--n_doped', type=int, default=1, help='Number of doping elements')
    parser.add_argument('-o', '--n_out', type=int, default=1, help='Number of output structures')
    parser.add_argument('-s', '--index_start', type=int, default=1, help='Starting index of the output structure filenames')
    parser.add_argument('-p', '--prefix', type=str, default='doped_structure_', help='Prefix of the output structure filenames')
    parser.add_argument('-l', '--existed_strus_dir_name', default='existed_structures', help='Directory name containing existing structures')

    args = parser.parse_args()

    input_file = args.input
    doped_element = args.doped_element
    doping_element = args.doping_element
    n_doped = args.n_doped
    n_out = args.n_out
    index_start = args.index_start
    prefix = args.prefix
    existed_strus_dir_name = args.existed_strus_dir_name
    
    if not os.path.exists(existed_strus_dir_name): os.makedirs(existed_strus_dir_name)
    
    tem_stru = Sr.from_file(input_file)  # Load the template structure from the input file 
    
    all_doping_strus = get_all_doping_strus(tem_stru, doped_element, doping_element, n_doped, n_out, index_start, prefix, existed_strus_dir_name)

    logging.info(f"Generated {len(all_doping_strus)} doped structures.")
    


if __name__ == '__main__':
    starttime = time.time()
    main()
    endtime = time.time()
    runtime = endtime-starttime
    logging.info("\nEnd of calculation.")
    logging.info("Program was running for %.2f seconds." % runtime)