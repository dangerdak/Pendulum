#!/bin/bash
./Pendulum
gnuplot -p << EOF
set datafile separator ','
set autoscale
set output "$1_energies.png"
set terminal png 

set format y ""
set ylabel "Energy"
set xlabel "Time"
set grid

plot "$1" using 1:4 title "dt = $4" with lines, "$2" using 1:4 title "dt = $5" with lines, "$3" using 1:4 title "dt = $6" with lines

reset

EOF
