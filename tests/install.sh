#!/bin/bash

# compile all files in lmpep_src
gfortran -c -O3 ../modules/rand_poly.f90
gfortran -c -O3 ../modules/initial_estimates.f90
gfortran -c -O3 ../modules/pzeros.f90

# create library
ar crv libmodules.a *.o
cp libmodules.a src/

# compile/run test_modules.f90
gfortran -O3 test_modules.f90 -L/Users/thcameron/Documents/FPML/tests/src -lmodules -o test_modules
./test_modules 1000 1000000 100

# compile test_modules.tex
pdflatex -shell-escape test_start_timing.tex
cp *.png figures/

#clean up
rm *.mod
rm *.o
rm *.a
rm *.log
rm *.aux
rm *.pdf
rm *.png
