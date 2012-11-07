#!/bin/bash
./Pendulum
gnuplot -p << EOF
set datafile separator ','
set title "Euler Energy Using Various Step Sizes"
set autoscale
set output "euler_energies.ps"
set terminal postscript

set ylabel "Energy"
set xlabel "time"
set grid

plot "$1" using 1:4 with lines, "$2" using 1:4 with lines,"$3" using 1:4 with lines

reset

EOF
