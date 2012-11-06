#!/bin/bash
./Pendulum
gnuplot -p << EOF
set datafile separator ','
set term x11 0
plot "$1" using 1:4 with lines, "$2" using 1:4 with lines,"$3" using 1:4 with lines
EOF
