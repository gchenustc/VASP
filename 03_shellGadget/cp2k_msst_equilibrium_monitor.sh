#!/bin/bash
# Author: gchen
# Time: 20241112
# Description: This script processes molecular dynamics simulation data of "CP2K" and uses gnuplot to generate plots for temperature over time, kinetic energy over time, and cell volume over frames. The script requires two input files: a data file (eg: project.ener) and a CP2K output file.
# Usage: ./your_script_name.sh <data_file> <cp2k_out_file>
# NOTE: If gnuplot fails to execute, the script will report an error and terminate.


# 检查参数是否存在
if [ $# -ne 2 ]; then
  echo "Usage: $0 <data_file> <cp2k_out_file>"
  exit 1
fi

# 获取参数
DATA_FILE=$1
CP2K_OUT_FILE=$2

# 检查数据文件是否存在
if [ ! -f "$DATA_FILE" ]; then
  echo "Error: Data file <DATA_FILE> not found!"
  exit 1
fi

if [ ! -f "$CP2K_OUT_FILE" ]; then
  echo "Error: Data file <CP2K_OUT_FILE> not found!"
  exit 1
fi

# 生成体积随步数变化的数据
grep "MD| Cell volume \[ang^3\]" $CP2K_OUT_FILE  | tr -s " " | cut -d " " -f6 | nl > volume-frame.txt

# 生成在x轴方向的压力随步数变化的数据
grep "P(xx) \[bar\]" $CP2K_OUT_FILE  | awk '{print $4}' | nl  > P_xx-frame.txt

# 生成在x轴方向的长度随步数变化的数据
grep "Cell lengths \[ang\]" $CP2K_OUT_FILE  | awk '{print $5}' | nl > x_lenth-frame.txt

# 生成 gnuplot 脚本文件 plot.gp
cat > plot.gp <<EOF

#############

# Temperature - Time
set terminal png enhanced font "Arial,12"
set output 'T.png'
set title "Temperature - Time"

set xlabel "Time (fs)"
set ylabel "Temperature (K)"

plot '$DATA_FILE' u 2:4 w l lw 2 t "Temperature - Time"

# Kinetic energy - Time
set terminal png enhanced font "Arial,12"
set output 'kin-E.png'
set title "Kinetic - Time"

set xlabel "Time (fs)"
set ylabel "Kinetic Energy (a.u.)"

plot '$DATA_FILE' u 2:3 w l lw 2 t "Kinetic - Time"

# Volume - Frame
set terminal png enhanced font "Arial,12"
set output 'vol.png'
set title "Volume - Frame"

set xlabel "Frames"
set ylabel "Volume (ang^3)"

plot 'volume-frame.txt' u 1:2 w l lw 2 t "Volume - Frame"


# P(xx) & X lenth - Frame
set terminal png enhanced font "Arial,12"
set output 'P-xx_lenth-x.png'
set y2tics

set title "P(xx) | X lenth - Frame"
set xlabel "Frames"
set ylabel "P(xx) (bar)"
set y2label "X lenth (Å)"

plot 'P_xx-frame.txt' u 1:2 w l lw 2 t "P(xx) - Frame" axis x1y1, \
     'x_lenth-frame.txt' u 1:2 w l lw 2  t "X lenth - Frame" axis x1y2

EOF

# 在shell中直接显示
cat > plot-dumb.gp <<EOF

# Temperature - Time
set terminal dumb
set title "Temperature - Time"

set xlabel "Time (fs)"
set ylabel "Temperature (K)"

plot '$DATA_FILE' u 2:4 w l t "Temperature - Time"

# Kinetic energy - Time
set terminal dumb
set title "Kinetic - Time"

set xlabel "Time (fs)"
set ylabel "Kinetic Energy (a.u.)"

plot '$DATA_FILE' u 2:3 w l t "Kinetic - Time"

# Volume - Frame
set terminal dumb
set title "Volume - Frame"

set xlabel "Frames"
set ylabel "Volume (ang^3)"

plot 'volume-frame.txt' u 1:2 w l t "Volume - Frame"

# P(xx) & lenth_x - Frame
set terminal dumb 
set y2tics

set title "P(xx) & X lenth - Frame"
set xlabel "Frames"
set ylabel "P(xx) (bar)"
set y2label "X lenth (Å)"

plot 'P_xx-frame.txt' u 1:2 w l t "P(xx) - Frame" axis x1y1, \
     'x_lenth-frame.txt' u 1:2 w l t "X lenth - Frame" axis x1y2

EOF

#############

# 调用 gnuplot 执行脚本
gnuplot plot.gp
gnuplot plot-dumb.gp

# 检查 gnuplot 是否成功执行
if [ $? -eq 0 ]; then
  echo "Plotting completed successfully. Output saved to output.png."
else
  echo "Error: Gnuplot execution failed."
  exit 1
fi

rm plot.gp
rm plot-dumb.gp
