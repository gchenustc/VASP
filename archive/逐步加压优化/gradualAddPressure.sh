# 设定压力
pressure="20 40 60 80 100"
current=-1
for pre in ${pressure}
do
    next=${pre}

    if [ -d ${next} ];then
        rm -f ${next}/*
    else
        mkdir ${next}
    fi
    cp INCAR POSCAR POTCAR KPOINTS sub.sh ./${next}

    if [ ${current} -ne -1 ];then
        cp ${current}/CONTCAR ${next}/POSCAR
    fi

    cd ${next}
    sed -i /PSTRESS/cPSTRESS=${next}0 INCAR
    srun vasp_std
    current=${next}
    cd ..
done

