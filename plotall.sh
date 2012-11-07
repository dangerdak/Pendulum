#!/bin/bash
gnuplot << EOF
set datafile separator ','
set autoscale
set nokey
set grid
set output "$1_plots.ps"
set terminal postscript

set size 1,1
set origin 0,0

set multiplot

set size 0.5,0.5
set xlabel "time"

#plot theta
set origin 0,0.5
set title "Theta over Time"
set ylabel "theta
plot "$1" u 1:2 with dots

#plot omega
set origin 0.5,0.5
set title "Omega over Time"
set ylabel "omega"
plot "$1" u 1:3 with dots

#plot energy
set origin 0,0
set title "Energy over Time"
set ylabel "Energy"
plot "$1" u 1:4 with dots

unset multiplot

reset

EOF
