#PBS -N N
#PBS -l nodes=1:ppn=28
#PBS -q paraque
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


for ((j=1; j<=50; j++)); do
    mkdir ${j}
    cp INCAR ${j}
    cp POTCAR ${j}
    cp submit.sh ${j}
    cp UCell_${j}_* ${j}/POSCAR

    cd ${j}
        mpirun -genv I_MPI_DEVICE ssm -machinefile /tmp/nodefile.$$ -n $NP /opt/issp2/vasp/vasp5.4.4_std
    cd ..

done

rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$
