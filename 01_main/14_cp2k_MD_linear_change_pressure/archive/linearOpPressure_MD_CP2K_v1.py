import os
import subprocess
import shutil

# 初始参数
initial_pressure = 100000  # bar - 1 GPa = 10000 bar
final_pressure = 500000  # bar
pressure_step = 50000
steps_per_segment = 1

total_steps = int((((final_pressure - initial_pressure) /  pressure_step ) + 1) * steps_per_segment)
project_name = "input"

# 初始结构
init_subsys = """
  &SUBSYS
    &CELL
      A    3.8     0.00000000     0.00000000
      B     0.00639460    3.8     0.00000000
      C     0.00000000     0.00000000    3.8
      PERIODIC XYZ #Direction(s) of applied PBC (geometry aspect)
    &END CELL
    &COORD
 N    0.324414927         0.324414955         0.324415239
 N    1.576953302         3.478321282         2.225783046
 N    3.478321595         2.225783272         1.576953135
 N    2.225783246         1.576953362         3.478321736
 N    2.225783272         2.225783272         2.225783272
 N    3.478321902         1.576953135         0.324415040
 N    1.576952602         0.324414955         3.478321736
 N    0.324415407         3.478321509         1.576953135
   &END COORD
#   &VELOCITY #You can set initial atomic velocities in this section
#   &END VELOCITY
    &KIND N
      ELEMENT N
      BASIS_SET DZVP-MOLOPT-SR-GTH-q5
      POTENTIAL GTH-PBE
    &END KIND
  &END SUBSYS
"""
  
template = """
&GLOBAL
  PROJECT {project}
  PRINT_LEVEL LOW
  RUN_TYPE MD
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
{subsys}
  &DFT
    BASIS_SET_FILE_NAME  BASIS_MOLOPT
    POTENTIAL_FILE_NAME  POTENTIAL
#   WFN_RESTART_FILE_NAME input-RESTART.wfn
    CHARGE    0 #Net charge
    MULTIPLICITY    1 #Spin multiplicity
    &QS
      EPS_DEFAULT 1.0E-10 #Set all EPS_xxx to values such that the energy will be correct up to this value
      EXTRAPOLATION ASPC #Extrapolation for wavefunction during e.g. MD. ASPC is default, PS can also be used
      EXTRAPOLATION_ORDER 3 #Order for PS or ASPC extrapolation. 3 is default
    &END QS
    &POISSON
      PERIODIC XYZ #Direction(s) of PBC for calculating electrostatics
      PSOLVER PERIODIC #The way to solve Poisson equation
    &END POISSON
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
#      &VDW_POTENTIAL
#        POTENTIAL_TYPE PAIR_POTENTIAL
#        &PAIR_POTENTIAL
#          PARAMETER_FILE_NAME dftd3.dat
#          TYPE DFTD3(BJ)
#          REFERENCE_FUNCTIONAL PBE
#          #CALCULATE_C9_TERM T #Calculate C9-related three-body term, more accurate for large system
#        &END PAIR_POTENTIAL
#      &END VDW_POTENTIAL
    &END XC
    &MGRID
      CUTOFF  300
      REL_CUTOFF  40
    &END MGRID
    &SCF
      MAX_SCF 25 #Maximum number of steps of inner SCF
      EPS_SCF 1.0E-05 #Convergence threshold of density matrix of inner SCF
#     SCF_GUESS RESTART #Use wavefunction from WFN_RESTART_FILE_NAME file as initial guess
#     IGNORE_CONVERGENCE_FAILURE #Continue calculation even if SCF not converged, works for version >= 2024.1
      &OT
        PRECONDITIONER FULL_KINETIC #FULL_SINGLE_INVERSE is also worth to try. FULL_ALL is better but quite expensive for large system
        MINIMIZER DIIS #CG is worth to consider in difficult cases
        LINESEARCH 2PNT #1D line search algorithm for CG. 2PNT is default. 3PNT is more expensive but may be better. GOLD is best but very expensive
        ALGORITHM STRICT #Algorithm of OT. Can be STRICT (default) or IRAC
      &END OT
      &OUTER_SCF
        MAX_SCF 20 #Maximum number of steps of outer SCF
        EPS_SCF 1.0E-05 #Convergence threshold of outer SCF
      &END OUTER_SCF
      &PRINT
        &RESTART OFF #Do not generate wfn file to suppress meaningless I/O cost
        &END RESTART
      &END PRINT
    &END SCF
#    SURFACE_DIPOLE_CORRECTION T
#    SURF_DIP_DIR Z
  &END DFT
  STRESS_TENSOR ANALYTICAL #Compute full stress tensor analytically
&END FORCE_EVAL

&MOTION
  &MD
    ENSEMBLE NPT_F
    STEPS {steps} #Number of steps to run
    TIMESTEP 1.0 #Step size in fs. Decrease it properly for high temperature simulation
    TEMPERATURE 300 #Initial and maintained temperature (K)
    COMVEL_TOL 0 #Uncomment this can remove translation motion of center-of-mass every step
    &THERMOSTAT
      TYPE CSVR
      &CSVR
        TIMECON  100 #Time constant in fs. Smaller/larger results in stronger/weaker temperature coupling
      &END CSVR
    &END THERMOSTAT
    &BAROSTAT
      PRESSURE {pressure} #Initial and maintained pressure (bar)
      TIMECON 1000 #Barostat time constant (fs)
      VIRIAL XYZ #Relax the cell along which cartesian axes
    &END BAROSTAT
    &PRINT
      &PROGRAM_RUN_INFO
        &EACH
          MD     1 #Output frequency of MD information, 0 means never
        &END EACH
      &END PROGRAM_RUN_INFO
    &END PRINT
  &END MD
#  &CONSTRAINT
#    &FIXED_ATOMS #Set atoms to be fixed
#      COMPONENTS_TO_FIX XYZ #Which fractional components will be fixed, can be X, Y, Z, XY, XZ, YZ, XYZ
#      LIST       1..320
#    &END FIXED_ATOMS
#  &END CONSTRAINT
  &PRINT
    &TRAJECTORY
      &EACH
        MD   1 #Output frequency of coordinates, 0 means never
      &END EACH
      FORMAT pdb
    &END TRAJECTORY
    &VELOCITIES
      &EACH
        MD     0 #Output frequency of velocities, 0 means never
      &END EACH
    &END VELOCITIES
    &FORCES
      &EACH
        MD     0 #Output frequency of forces, 0 means never
      &END EACH
    &END FORCES
    &RESTART
      BACKUP_COPIES 0 #Maximum number of backing up restart file, 0 means never
      &EACH
        MD  1 #Frequency of updating last restart file, 0 means never
      &END EACH
    &END RESTART
    &RESTART_HISTORY OFF 
    &END RESTART_HISTORY
  &END PRINT
&END MOTION
"""

def extract_subsys(path):
    # 定义文件路径
    file_path = path

    # 读取文件内容
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # 初始化标志变量
    in_subsys = False
    extracted_lines = []

    # 遍历每一行
    for line in lines:
        if line.strip().startswith('&SUBSYS'):
            in_subsys = True
            extracted_lines.append(line)
        elif line.strip().startswith('&END SUBSYS'):
            in_subsys = False
            extracted_lines.append(line)
        elif in_subsys:
            extracted_lines.append(line)

    # 合并提取的行
    extracted_text = ''.join(extracted_lines)

    if extracted_text:
        return extracted_text
    else:
        print('No matching text found.')

# 初始压力
current_pressure = initial_pressure

for step in range(0, total_steps, steps_per_segment):

    # 更新输入文件
    if current_pressure == initial_pressure:
        input_content = template.format(steps=steps_per_segment, pressure=current_pressure, project=project_name, subsys=init_subsys)
    else: # 当非初始压力时，寻找上一步计算完成的结构段落 &subsys
        last_pressure = current_pressure - pressure_step
        last_pressure_path = f"{int(last_pressure)}_bar/{project_name}-1.restart"
        current_subsys = extract_subsys(last_pressure_path)
        input_content = template.format(steps=steps_per_segment, pressure=current_pressure, project=project_name, subsys=current_subsys)
        
    # 创建当前压力下的文件夹
    os.makedirs(f"{int(current_pressure)}_bar", exist_ok=True)
    os.chdir(f"{int(current_pressure)}_bar")
    
    # 创建inp
    with open(f'{project_name}.inp', 'w') as f:
        f.write(input_content)
    # 创建其他相关文件
    shutil.copy("/public/home/c1337375425/storage/cp2k/BASIS_MOLOPT", "./")
    shutil.copy("/public/home/c1337375425/storage/cp2k/POTENTIAL", "./")
    
    # 运行 CP2K
    subprocess.run(['srun', "--mpi=pmi2", "cp2k.psmp", "-i", f'{project_name}.inp'])
    
    # 退出当前文件夹
    os.chdir(f"../")
    
    # 更新压力
    current_pressure += pressure_step

print("Simulation completed.")
