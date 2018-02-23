#!/bin/bash

# compile all files in lmpep_src
gfortran -c ../modules/rand_poly.f90
gfortran -c ../modules/initial_estimates.f90
gfortran -c ../modules/modified_laguerre.f90
gfortran -c ../modules/pzeros.f90
gfortran -c ../modules/fpml.f90

# create library
ar crv libmodules.a *.o
cp libmodules.a src/

# compile test_modules.f90
gfortran test_modules.f90 -L/Users/thcameron/Documents/FPML/tests/src -lmodules -o test_modules
gfortran test_main.f90 -L/Users/thcameron/Documents/FPML/tests/src -lmodules -o test_main

#clean up
rm *.mod
rm *.o
rm *.a
