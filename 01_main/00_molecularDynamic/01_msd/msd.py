"""
Desc:
    This script calculates the mean square displacement (MSD) of atoms in a molecular dynamics simulation and generates plots of the MSD. The MSD measures the average displacement of atoms from their initial positions over time.
Usage:
    python file.py -p 0.5 # -p is the length for each step
Author: 
    gchen
Time: 
    2022-11-4
"""

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib as mpl
import time


def plot():
    # Configuration
    mpl.use('TkAgg')
    plt.rcParams['font.sans-serif']=['Arial']
    formatter = ticker.FormatStrFormatter("%.2e")
    msd = pd.read_csv("msd.out",sep=" ",index_col="time(fs)")
    # Plot figure 1
    fig = plt.figure(dpi=60,figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.plot(msd.index,msd.iloc[:,0:4],aa=True)
    ax.legend(msd.columns,frameon=False,loc="upper center",ncol=4)
    _extracted_from_plot_12(ax, msd, formatter, "msd.png")
    # Plot figure 2
    fig1,ax1 = plt.subplots(1,1,dpi=60,figsize=(10,8))
    ax1.plot(msd.index,msd.iloc[:,4:],aa=True)
    ax1.legend(msd.iloc[:,4:].columns,frameon=False,loc="upper center",ncol=len(msd.iloc[:,4:].columns))
    _extracted_from_plot_12(ax1, msd, formatter, "msd-element.png")


def _extracted_from_plot_12(arg0, msd, formatter, arg3):
    arg0.set_title("MSD", fontsize=15)
    arg0.set_xlabel(msd.index.name)
    arg0.set_ylabel("MSD(A^2)")
    arg0.set_xlim(0, msd.index[-1])
    arg0.yaxis.set_major_formatter(formatter)
    plt.savefig(arg3)
    plt.close()

def car2fra(pos_arr, const):
    """Convert Cartesian coordinates to fractional coordinates.

    Args:
        pos_arr (np.ndarray): Cartesian coordinates.
        const (np.ndarray): Lattice constants.

    Returns:
        np.ndarray: Fractional coordinates.
    """
    return np.dot(pos_arr, np.linalg.inv(const))

def fra2car(pos_arr, const):
    """Convert fractional coordinates to Cartesian coordinates.

    Args:
        pos_arr (np.ndarray): Fractional coordinates.
        const (np.ndarray): Lattice constants.

    Returns:
        np.ndarray: Cartesian coordinates.
    """
    return np.dot(pos_arr,const)
 
def judge_period_atoms(before_pos,after_pos,cons,record):
    """Check the periodicity of atomic positions.

    Args:
        before_pos (np.ndarray): Positions before the update.
        after_pos (np.ndarray): Positions after the update.
        cons (np.ndarray): Lattice constants.
        record (np.ndarray): Record of atom positions beyond the boundary.
    """
    for i in range(3):
        pos_diff = after_pos[:,i] - before_pos[:,i]
        record[:,i] += np.where(pos_diff < -0.5,1,0)
        record[:,i] += np.where(pos_diff > 0.5,-1,0)
    
    
def calibration_pos(pos,record):
    """Calibrate positions based on periodic boundary conditions.

    Args:
        pos (np.ndarray): Positions.
        record (np.ndarray): Record of atom positions beyond the boundary.

    Returns:
        np.ndarray: Calibrated positions.
    """
    for i in range(3): #3d
        pos[:,i] = pos[:,i] + record[:,i]
    return pos
    
def calc_msd(file1,file2,demision=3,pmsd_list=None):
    """Calculate mean square displacement (MSD).


    Args:
        file1 (np.ndarray): Positions of atoms.
        file2 (np.ndarray): Reference positions.
        demision (int): Dimensionality of the calculation (1 or 3).
        pmsd_list (List[int]): List of atomic counts for partial MSD calculation.

    Returns:
        np.ndarray or List[np.ndarray]: MSD or partial MSD values.
    """
    assert file1.shape == file2.shape

    if pmsd_list:
        return [
            calc_msd(i, j)
            for i, j in zip(
                np.vsplit(file1, np.cumsum(pmsd_list)[:-1]),
                np.vsplit(file2, np.cumsum(pmsd_list)[:-1]),
            )
        ]
    if demision == 3:
        return np.square(file1-file2).sum(axis=1).mean()
    elif demision == 1:
        return np.square(file1-file2).mean(axis=0)

if __name__ == "__main__":
    start_time = time.perf_counter()

    parser = argparse.ArgumentParser(description='give step length')
    parser.add_argument("-p","--potim", action="store", type=float, required=True, help="give potim for incar")
    options = parser.parse_args()

    from ase.io.vasp import read_vasp_xdatcar
    xdatcar = read_vasp_xdatcar("XDATCAR",index=0)
    xdatcar_0 = xdatcar[0]

    # Atom types -- e.g., ['H', 'N']
    species = list(set(xdatcar_0.get_chemical_symbols()))
    # Number of atoms for each atom type -- e.g., [12, 192]
    n_atoms = [len([atom for atom in xdatcar_0 if atom.symbol==specie]) for specie in species]
    n_atoms_total = sum(n_atoms)
    # Step length
    optim = options.potim
    # Coordinates of the first structure in Cartesian coordinates
    init_pos = xdatcar_0.get_positions()
    # Dictionary of atom information, e.g., {"N":192,"H":8}
    atoms_info_dict={}
    for index,specie in enumerate(species):
        atoms_info_dict.setdefault(specie, n_atoms[index])

    # Output
    msd_list = ['time(fs)','MSD-total','MSD-x','MSD-y','MSD-z']  # Header line for MSD output file
    msd_list.extend('MSD-'+i for i in atoms_info_dict)
    # Record the number of times atoms exceed the boundary, the y-axis corresponds to each atom, and the x-axis corresponds to the x, y, z coordinates of each atom.
    over_boundary_record = np.zeros((n_atoms_total,3)) 

    pos_bak=None # Positions in the previous step
    for n_step,atoms in enumerate(xdatcar):
        cons = atoms.cell  # Lattice constants for each step
        pos = car2fra(atoms.get_positions(), cons)

        if n_step != 0:
            # Positions in the previous step
            judge_period_atoms(pos_bak, pos, cons, over_boundary_record)

        pos_bak = pos.copy()

        # Calibrate positions
        pos = calibration_pos(pos,over_boundary_record)
        pos_car = fra2car(pos, cons)

        msd_list_append = [
            (n_step + 1) * optim,
            calc_msd(pos_car, xdatcar_0.get_positions(), 3),
        ]
        msd_list_append.extend(calc_msd(pos_car,xdatcar_0.get_positions(),1)) # Partial MSD
        msd_list_append.extend(calc_msd(pos_car,xdatcar_0.get_positions(),3,list(atoms_info_dict.values()))) # Partial MSD
        msd_list  = np.vstack([msd_list,msd_list_append])

    np.savetxt("msd.out",msd_list,fmt="%s")
    plot()  # Plotting

    time_cost = time.perf_counter () - start_time
    print(
        f'---- Time used: {time.strftime("%H:%M:%S", time.gmtime(time_cost))} ----'
    )
