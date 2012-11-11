#!/bin/bash
./Pendulum
gnuplot -p << EOF
set datafile separator ','
set autoscale
set output "euler_energies.png"
set terminal png enhanced 

set ylabel "Energy"
set xlabel "Time"
set grid

plot "$1" using 1:4 title "dt = 0.199" with lines, "$2" using 1:4 title "dt = 0.200" with lines, "$3" using 1:4 title "dt = 0.201" with lines

reset

EOF
