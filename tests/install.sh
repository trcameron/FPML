#!/bin/bash

# compile all files in lmpep_src
gfortran -c -g -Wall ../modules/rand_poly.f90
gfortran -c -g -Wall ../modules/initial_estimates.f90
gfortran -c -g -Wall ../modules/modified_laguerre.f90
gfortran -c -g -Wall ../modules/pzeros.f90

# create library
ar crv libmodules.a *.o
cp libmodules.a src/

# compile test_modules.f90
gfortran -g -Wall test_modules.f90 -L/Users/thcameron/Documents/FPML/tests/src -lmodules -o test_modules

#clean up
rm *.mod
rm *.o
rm *.a
