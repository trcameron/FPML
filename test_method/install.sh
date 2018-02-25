#!/bin/bash

# compile all files in modules
gfortran -c -O2 modules/rand_poly.f90
gfortran -c -O2 modules/initial_estimates.f90
gfortran -c -O2 modules/methods.f90
gfortran -c -O2 modules/driver.f90

# create library
ar crv libtest_method.a *.o
cp libtest_method.a library/

# compile test_methods.f90
gfortran -O2 tests/test_methods.f90 -L/Users/thcameron/Documents/FPML/test_method/library -ltest_method -o test_method

# clean up
rm *.mod
rm *.o
rm *.a