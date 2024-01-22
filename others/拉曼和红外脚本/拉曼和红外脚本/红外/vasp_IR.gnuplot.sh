#!/bin/bash

gnuplot << EOF
reset
set terminal postscript enhanced "Helvetica" 20
set output "Raman.eps"
set xrange[0:4000]
set style line 1 lt 1 lc -1 lw 4 
set title "Off-resonance Raman PBE/PAW_h level"
#set nokey
set xlabel "Frequency, cm^{-1}"
set ylabel "Intensity"

plot "results.txt-broaden.dat" u 1:2 w l ls 1
EOF

ps2pdf IR.eps
