#!/apps/python3/3.10.0/bin/python3

"""
This script analyzes a CP2K output file and extracts information about the
geometry optimization, including the number of SCF cycles, energy change, and
convergence criteria.

Usage:
  cp2k-output-analyse.py OUTPUT_FILE.out
"""

import os,sys,datetime

os.chdir(os.getcwd())

print('''
***************************************
*** The format of the script: ***
cp2k-output-analyse.py OUTPUT_FILE.out
***************************************
      ''')

try:
    input_file = sys.argv[1]
except IndexError:
    print("\033[31mERROR:\033[0m Missing file operand! Please identify the name of OUTPUT_FILE.out")
    sys.exit() 

print("=======================================")
print("In process...")
print("...")


output_file = sys.argv[1]
plot_file = output_file.split('.out')[0]  + "__data.csv"

with open(output_file, 'r') as f:
    with open(plot_file,'w') as o:
        lines = f.readlines()
        num_lines = len(lines) # the number of lines in the file
        MAX_SCF = 50
        OUTER_SCF_CHECK = "FALSE"
        SCF_OPTIMIZER = "DIAGONALIZATION"
        GEO_OPTIMIZER = "N/A"
        TOTAL_TIME = 0
        for line in lines:
            if "PROGRAM STARTED AT" in line: 
                starttime = line.split('AT')[-1].strip() #use strip to delete the space in front and behind the string
            elif "Run type" in line:
                RUN_TYPE = line.split()[-1]
                if RUN_TYPE == "GEO_OPT":
                    line_added = 25
                elif RUN_TYPE == "CELL_OPT":
                    line_added = 32
                else:
                    print("\033[31mERROR:\033[0m This script can only be used for Geometry Optimization results!")
                    sys.exit()
            elif "eps_scf:" in line:
                EPS_SCF = line.split(':')[-1].strip()
            elif "Outer loop SCF in use" in line:
                OUTER_SCF_CHECK = "TRUE"
            elif "max_scf:" in line:
                MAX_SCF = line.split(':')[-1].strip()
            elif "STARTING GEOMETRY OPTIMIZATION" in line:
                line_number = lines.index(line)
                line_number += 1  # grab the next line
                GEO_OPTIMIZER = lines[line_number].split('***')[1].strip()
            elif " OT " in line:
                SCF_OPTIMIZER = "OT"
            elif "SCF run NOT converged ***" in line:
                SCF_NUMBER = "!" + MAX_SCF
            elif "*** SCF run converged" in line:
                SCF_NUMBER = line.split()[-3]
            elif "outer SCF loop converged in" in line:
                SCF_NUMBER = line.split()[-2]
            elif "Informations at step" in line:
                CYCLE_NUMBER = line.split('=')[-1].split('-')[0].strip()
                top_line_number = lines.index(line)
                bottom_line_number = top_line_number + line_added 
                for contents in lines[top_line_number:bottom_line_number]:
                    if "Decrease in energy " in contents:
                        ENERGY_CHANGE = contents.split('=')[-1].strip() 
                    elif "Total Energy" in contents:
                        TOTAL_ENERGY = contents.split('=')[-1].strip()
                    elif "Conv. limit for step size" in contents:
                        MAX_D = "0" + contents.split('=')[-1].strip().strip('0')
                    elif "Conv. limit for RMS step" in contents:
                        RMS_D = "0" + contents.split('=')[-1].strip().strip('0')
                    elif "Conv. limit for gradients" in contents:
                        MAX_F = "0" + contents.split('=')[-1].strip().strip('0')
                    elif "Conv. limit for RMS grad" in contents:
                        RMS_F = "0" + contents.split('=')[-1].strip().strip('0')
                    elif "Max. step size " in contents:
                        MAX_D_VALUE = float(contents.split('=')[-1].strip())
                    elif "RMS step size " in contents:
                        RMS_D_VALUE = float(contents.split('=')[-1].strip())
                    elif "Max. gradient " in contents:
                        MAX_F_VALUE = float(contents.split('=')[-1].strip())
                    elif "RMS gradient " in contents:
                        RMS_F_VALUE = float(contents.split('=')[-1].strip())
                    elif "Used time" in contents:
                        USEDTIME = contents.split('=')[-1].strip()
                        TOTAL_TIME += float(USEDTIME)
                try:
                    ENERGY_CHANGE = ENERGY_CHANGE
                    if ENERGY_CHANGE == "NO":
                        o.write("%2s %4s %19s %7s"%("xx",CYCLE_NUMBER,TOTAL_ENERGY,SCF_NUMBER))
                    else:
                        o.write("%7s %19s %7s"%(CYCLE_NUMBER,TOTAL_ENERGY,SCF_NUMBER))
                except NameError:
                    o.write("%7s %19s %7s"%(CYCLE_NUMBER,TOTAL_ENERGY,SCF_NUMBER))
                try:
                    o.write("%12.6f"%(MAX_D_VALUE))
                    if MAX_D_VALUE > float(MAX_D):
                        MAX_D_CONVERGENCE = "NO"
                    else:
                        MAX_D_CONVERGENCE = "YES"
                    o.write("%4s"%(MAX_D_CONVERGENCE))
                    o.write("%12.6f"%(RMS_D_VALUE))
                    if RMS_D_VALUE > float(RMS_D):
                        RMS_D_CONVERGENCE = "NO"
                    else:
                        RMS_D_CONVERGENCE = "YES"
                    o.write("%4s"%(RMS_D_CONVERGENCE))
                    o.write("%14.6f"%(MAX_F_VALUE))
                    if MAX_F_VALUE > float(MAX_F):
                        MAX_F_CONVERGENCE = "NO"
                    else:
                        MAX_F_CONVERGENCE = "YES"
                    o.write("%4s"%(MAX_F_CONVERGENCE))
                    o.write("%13.6f"%(RMS_F_VALUE))
                    if RMS_F_VALUE > float(RMS_F):
                        RMS_F_CONVERGENCE = "NO"
                    else:
                        RMS_F_CONVERGENCE = "YES"
                    o.write("%4s"%(RMS_F_CONVERGENCE))
                    o.write("%13s"%(USEDTIME))
                except NameError:
                    o.write("%11s %15s %16s %17s"%("N/A", "N/A", "N/A", "N/A"))
                    o.write("%18s"%(USEDTIME))
                o.write("\n")
            elif "OPTIMIZATION COMPLETED" in line:
                CYCLE_NUMBER = line.split('=')[-1].split('-')[0].strip()
                top_line_number = lines.index(line)
                bottom_line_number = num_lines
                for contents in lines[top_line_number:bottom_line_number]:
                    if "ENERGY" in contents:
                        TOTAL_ENERGY = contents.split(':')[-1].strip()
                        TOTAL_ENERGY = float(TOTAL_ENERGY)
                        TOTAL_ENERGY = "%.10f" % TOTAL_ENERGY
                
                o.write("%9s %17s %7s"%("Final",str(TOTAL_ENERGY),SCF_NUMBER))
                o.write("%11s %15s %16s %17s"%("N/A", "N/A", "N/A", "N/A"))
                o.write("%17s"%("N/A"))
                o.write("\n")
        TOTAL_TIME = str(datetime.timedelta(seconds=float(TOTAL_TIME)))
        o.write("# Done!")
        o.close()
    f.close()

with open(plot_file, 'r+') as f:
    contents = f.read()
    f.seek(0, 0)
    f.write("# Job Starting Date: " + starttime \
            + "\n# Total used time: " + str(TOTAL_TIME) \
            + "\n# Directory: " + os.getcwd() \
            + "\n# RUN_TYPE: " + RUN_TYPE \
            + "\n# EPS_SCF: " + EPS_SCF \
            + "\n# MAX_SCF: " + MAX_SCF \
            + "\n# SCF_OPTIMIZER: " + SCF_OPTIMIZER \
            + "\n# OUTER_SCF: " + OUTER_SCF_CHECK \
            + "\n# GEO_OPTIMIZER: " + GEO_OPTIMIZER \
            + "\n# CYCLE | TOTAL_ENERGY [a.u.] | SCF | MAX_D.(" + MAX_D + ") | RMS_D.(" + RMS_D + ") | MAX_F.(" + MAX_F + ") | RMS_F.(" + RMS_F + ") | USEDTIME [s]" \
            + "\n" + contents)
    f.close()
print("=======================================")


print("Do you want to plot cycle vs. energy?\n(y/n)")
plot_choice = input()
if plot_choice == "y":
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import MultipleLocator

    f = open(plot_file, 'r')
    x = []
    y = []
    for lines in f:
        x_value = lines.split()[0]
        if x_value.isdigit():
            x.append(int(x_value))
            y.append(float(lines.split()[1]))
        elif x_value == "xx":
            x.append(int(lines.split()[1]))
            y.append(float(lines.split()[2]))

    plt.scatter(x, y)
    plt.xlabel("Cycle Number")
    plt.ylabel("Energy (a.u.)")
    if len(x) <= 50:
        x_spacing = 5
    elif len(x) <= 100:
        x_spacing = 10
    elif len(x) <= 150:
        x_spacing = 15
    elif len(x) <= 200:
        x_spacing = 20
    else:
        x_spacing = 50
    x_major_locator=MultipleLocator(x_spacing)
    ax = plt.gca()
    ax.xaxis.set_major_locator(x_major_locator)
    plt.show()
print("Done!")
print("=======================================")






