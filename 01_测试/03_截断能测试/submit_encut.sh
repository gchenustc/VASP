#!/bin/env bash
#PBS -N 100a
#PBS -l nodes=1:ppn=40
#PBS -q paraque3
#PBS -V
#PBS -o ./
#PBS -e ./
#PBS -S /bin/bash 
#### Set intel environment###
source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/intel64/mklvars_intel64.sh
source /opt/intel/impi/5.0.2.044/bin64/mpivars.sh

cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

#############需要修改的参数###################
list="300 320 340 360 380 400 420 440 460 480 500 520 540 560 580 600 620 640 660 680 700 720 740 760 780 800"
#############需要修改的参数###################

for i in $list
do
    workDir=Encut_$i
    [ -d $workDir ] || echo "the target dir not existexit"\>$workDir > error
    cd $workDir
    mpirun -genv I_MPI_DEVICE ssm -machinefile /tmp/nodefile.$$ -n $NP /opt/issp2/vasp/vasp5.4.4_std
    cd ..
done

rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$
