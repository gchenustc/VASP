# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 00:24:14 2018, version 3
@author: jliubt
Contact: suntan1991@gmail.com
"""

# This script is inpsired by previous one developed by Muratahan Aykol
# To accurately calculate the msd with block sampling
# Every neighbouring time steps should check with nearest images considering pbc
# It is better to use this code for total time steps larger than 10,000. This is
# very common in any MD simulations, both classical and ab-initio MD
import numpy as np
import os, pickle, time


def MSD(xdata_file):
    
    xdatcar = open(xdata_file, 'r')
    lines = xdatcar.readlines()
    xdatcar.close()
    len_xdatcar = len(lines)                                                    # get the total lenght of the XDATCAR file to decide the time step of MD
    print (len_xdatcar)
    
    xdatcar = open(xdata_file, 'r')
    system = xdatcar.readline()
    print ("Calculating Mean Squared Displacement (MSD) for system %s" % system)
    
    # get the scale and the lattice vectors in the head of the XDATCAR file
    scale = float(xdatcar.readline().rstrip('\n'))
    lattice = np.zeros((3,3))
    #basis vectors in cartesian coords
    for i in range(3):
        lattice[i] = np.array([float(s)*scale for s in xdatcar.readline().rstrip('\n').split()])

    # get the element of the system and the correpsonding atom numbers
    element_list = xdatcar.readline().rstrip('\n').split()
    element_dict = {}
    element_numbers = [int(s) for s in xdatcar.readline().rstrip('\n').split()]
    N_atoms = sum(element_numbers)
    
    for t in range(len(element_list)):
        element_dict[element_list[t]] = element_numbers[t]
    
    # MD time steps calculated
    N_frame = (len_xdatcar-7) // (N_atoms+1)
    print ("MD simulation runs for total %d time steps" % N_frame)
    
    # Store all the positions in the variable pos
    if os.path.isfile('xdatcar.pckl'):
        f = open('xdatcar.pckl', 'rb')
        pos = pickle.load(f)
        f.close()
    else:
        pos = np.zeros((N_frame, N_atoms, 3))
        for t in range(N_frame):
            xdatcar.readline()
            for j in range(N_atoms):
                pos[t,j] = np.array([ float(s) for s in xdatcar.readline().rstrip('\n').split() ])
            if t % 1000 == 0:
                print (t)
        f = open("xdatcar.pckl", "wb")
        pickle.dump(pos, f)
        f.close()
    xdatcar.close()                                                             # close the file
    print ("Successfully read the XDATACAR file into the program!!!")
    
    recorder = open("msd_jp.out", 'w')
    recorder.write("step ")
    for element in element_list:
        recorder.write(element+" ")
    recorder.write("\n")
    
    ###########################################################################
    # going to calculate the Mean-Squared-Displacement
    # minimum unit steps for msd sampling, this can be adjusted if needed
    block = 1
    N_block = (N_frame-1) // block
    msd_block = np.zeros((N_block, N_atoms))
    msd_lattice_block = np.zeros((N_block, N_atoms,3))
    
    for k, step in enumerate(range(block, N_frame, block)):
        msd_sample = np.zeros((N_atoms))
        msd_lattice_sample = np.zeros((N_atoms,3))
        
        # divide the total time steps into small blocks, one block with [step] images
        # and sum all the images in subblock, then average
        # Here I use two differnt schemes for sampling: 
        # get the distance btw two cascade time steps, and then add N_block distances together to get msd
        count = 0
        for j in range(0, step, 10):
            pos_copy = []
            pos_copy = pos[j::step]                                             # Sampling for every [step]th entry 
            diff_pos = np.diff(pos_copy, axis=0)                                # distance in fractional coordinate
            #get nearest distance consiering pbc, e.g., 0.2-->0.2; 0.8-->-0.2; -0.3-->-0.3; -0.7-->0.3
            dist = diff_pos - np.rint(diff_pos)
            dist = np.dot(dist,lattice)                                         # distance in cartesian coordinate
            dist_square = dist*dist
            msd_lattice_sample += np.sum(dist_square,axis=0)                    # take average of different sampling
            msd_sample += np.sum(np.sum(dist_square,axis=0), axis=1)
            count += len(diff_pos)
        msd_lattice_sample /= count
        msd_sample /= count

        if step % 100 == 0:
            print (step)
        msd_block[k] = msd_sample
        msd_lattice_block[k] = msd_lattice_sample
        
        
        recorder.write(str(step) + " ")
        el_index_start = 0
        for el in element_list:
            el_index_end = el_index_start + element_dict[el]
            msd_el = np.mean(msd_block[k][el_index_start:el_index_end])
            msd_lattice_el = np.mean(msd_lattice_block[k][el_index_start:el_index_end], axis=0)
            recorder.write(str(msd_el)+ " " + str(msd_lattice_el[0])+ " " + str(msd_lattice_el[1])+ " " + str(msd_lattice_el[2])+ " ")
            el_index_start = el_index_end
        recorder.write("\n")
    recorder.close()

start_time = time.perf_counter()
MSD("XDATCAR")
time_cost = time.perf_counter () - start_time
print("------ Time used: %s ------" % time.strftime("%H:%M:%S", time.gmtime(time_cost)))