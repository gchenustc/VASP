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
list=("1 1 1" "2 2 2" "3 3 3" "4 4 4" "5 5 5" "6 6 6" "7 7 7" "8 8 8" "9 9 9" "10 10 10" "11 11 11" "12 12 12" "13 13 13")
#############需要修改的参数###################

KpCount=${#list[*]}  # 统计测试的k点个数
for i in $(seq 0 1 $((${KpCount}-1))) 
do
    workDir=kpoints_`echo $(echo ${list[$i]}) | awk 'BEGIN{OFS="-"};{print $1,$2,$3}'`
    [ -d $workDir ] || echo "the target dir not existexit"$workDir > error
    cd $workDir
    mpirun -genv I_MPI_DEVICE ssm -machinefile /tmp/nodefile.$$ -n $NP /opt/issp2/vasp/vasp5.4.4_std
    cd ..
done

rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$
