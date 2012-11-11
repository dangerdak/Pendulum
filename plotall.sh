#!/bin/bash
gnuplot << EOF
set datafile separator ','
set autoscale
set grid
set output "$1_plots.png"
set terminal png font "sans, 7"

set size 1,1
set origin 0,0

set multiplot

set size 0.5,0.5
set xlabel "time"
set ylabel "Energy"

#plot theta
set origin 0,0.5
plot "$1" u 1:4 title "$4" with lines

#plot omega
set origin 0.5,0.5
plot "$2" u 1:4 title "$5" with lines

#plot energy
set origin 0,0
plot "$3" u 1:4 title "$6" with lines

unset multiplot

reset

EOF
