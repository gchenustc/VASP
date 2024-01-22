#!/bin/sh
#Desc: print the forces and energy on the screen. The first parameter is a number that determines how many steps are printed, and defaults to 1

begin=$1
if [ $begin = '' ]; then begin=1; fi    
awk -v begin=$begin '/E0/{if ( i<begin  ) i++;else print $0 }' OSZICAR >temp.e
    
awk '/POSITION/,/drift/{
    if(NF==6) print $4,$5,$6;
    else if($1=="total") print $1 }' OUTCAR >temp.f

awk '{if($4=="F"||$4=="T") print $4,$5,$6}' CONTCAR >temp.fix
flag=`wc temp.fix|awk '{print $1}'`
steps=`grep E0 OSZICAR |tail -1 |awk '{print $1}'`
if [ flag != '0' ] ; then
    if [ -f temp.fixx ] ; then rm temp.fixx ; fi
    for i in `seq $steps`;do
        cat temp.fix >>temp.fixx
        echo >>temp.fixx
    done
    paste temp.f temp.fixx >temp.ff
fi

awk  '{ if($1=="total") {print ++i,a;a=0}
        else {
            if($4=="F") x=0; else x=$1;
            if($5=="F") y=0; else y=$2;
            if($6=="F") z=0; else z=$3;
            force=sqrt(x^2+y^2+z^2);
            if(a<force) a=force} }' temp.ff >force.conv

echo  steps  force:
tail  -${1}  force.conv
echo
echo  steps  energy:
tail  -${1}  temp.e

rm temp.e temp.f temp.ff temp.fix temp.fixx force.conv
