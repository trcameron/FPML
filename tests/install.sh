#!/bin/bash

# compile all files in lmpep_src
gfortran -c -O2 ../modules/rand_poly.f90
gfortran -c -O2 ../modules/initial_estimates.f90
gfortran -c -O2 ../modules/modified_laguerre.f90
gfortran -c -O2 ../modules/pzeros.f90
gfortran -c -O2 ../modules/fpml.f90

# create library
ar crv libmodules.a *.o
cp libmodules.a src/

# compile test_modules.f90
gfortran -O2 test_modules.f90 -L/Users/thcameron/Documents/FPML/tests/src -lmodules -o test_modules
gfortran -O2 test_main.f90 -L/Users/thcameron/Documents/FPML/tests/src -lmodules -o test_main

#clean up
rm *.mod
rm *.o
rm *.a
