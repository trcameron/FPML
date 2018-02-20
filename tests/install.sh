#!/bin/bash

# compile all files in lmpep_src
gfortran -c -O3 ../modules/*

# create library
ar crv libmodules.a *.o
cp libmodules.a src/

# compile/run test_modules.f90
gfortran -O3 test_modules.f90 -L/Users/thcameron/Documents/FPML/tests/src -lmodules -o test_modules
./test_modules

# compile test_modules.tex
pdflatex -shell-escape test_modules.tex
cp *.png figures/

#clean up
rm *.mod
rm *.o
rm *.a
rm *.log
rm *.aux
rm *.pdf
rm *.png
