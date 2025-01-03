import dpdata
import os
import argparse
"""
DESC: convert the deepmd's data set into POSCARs.
USAGE:
        python /path/to/script.py -d dp_path -p poscar_stored_path
"""

def createDir(dirName):
    folder_name = dirName  

    if not os.path.exists(folder_name):  
        os.makedirs(folder_name)  


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Specify the destination folder path.")  
    parser.add_argument("-d", "--dp_path", type=str, help="Specify the deepmd data set folder path.", default="deepmd")  
    parser.add_argument("-p", "--poscar_path", type=str, help="Specify the folder path for store poscars.", default="POSCARs")  
    args = parser.parse_args()  

    POSCAR_PATH = args.poscar_path
    createDir(POSCAR_PATH)

    d_dp = dpdata.LabeledSystem(args.dp_path, fmt = 'deepmd/npy')
    for i in range(len(d_dp)):
        d_dp.to('vasp/poscar', './' + str(POSCAR_PATH) + '/POSCAR.' + str(i+1) + '.vasp' , frame_idx=i)