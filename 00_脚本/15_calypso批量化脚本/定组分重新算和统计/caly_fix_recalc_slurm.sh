for ((j=1; j<=50; j++)); do
    mkdir ${j}
    cp INCAR ${j}
    cp POTCAR ${j}
    cp submit.sh ${j}
    cp UCell_${j}_* ${j}/POSCAR

    cd ${j}
        sbatch submit.sh
    cd ..

done

