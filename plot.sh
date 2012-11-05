#!/bin/bash
./Pendulum
gnuplot -p << EOF
set datafile separator ','
set term x11 0
plot "$1" using 1:2 with lines
set term x11 1
plot "$1" using 1:3 with lines
set term x11 2
plot "$1" using 1:4 with lines
set term x11 3
# splot "$1" using 1:2:3 with lines
EOF
