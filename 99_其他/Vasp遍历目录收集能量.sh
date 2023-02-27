#!/bin/env bash
# desc: print enthalpy value for every pressure, such as the down level directory name are 100 200 300,respectively.  present 100,200,300GPa for vasp's structure relaxation, this file is in upper level directory."
# author: gchen
echo "time: `date`" > enthalpy.txt

for i in {140..290..10} ;do
cd ${i}
echo $(awk -F" " -v label=${i} 'BEGIN{OFS="\t"};{if($4=="E0="){print label,$5}}' OSZICAR | tail -1) >> ~/Task/H-P-N/BlueN/enthalpy.txt
cd ..
done
