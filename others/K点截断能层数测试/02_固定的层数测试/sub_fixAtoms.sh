#PBS -N fix
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

fix_list="6 12"  #固定的测试层数
for fix_num in $fix_list
do
    # 创建文件夹
    workDir=fix-${fix_num}atoms
    [ -d ${workDir} ] || echo "the destinate dir not exist"

    cd ${workDir}

    mpirun -genv I_MPI_DEVICE ssm -machinefile /tmp/nodefile.$$ -n $NP /opt/issp2/vasp/vasp5.4.4_std

    cd ..
done

rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$
