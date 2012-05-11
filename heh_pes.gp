#!/usr/bin/env gnuplot
set xzeroaxis lt -1
set xlabel "Distance (Bohr)"
set ylabel "Energy (Hartree)"
plot "heh_pes.dat" smooth csplines
