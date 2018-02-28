#!/bin/bash

# compile AMVW_src
gfortran -c -O2 AMVW_src/*.f90
# create library
ar crv libAMVW.a *.o
cp libAMVW.a AMVW_src/

# compile test_methods.f90
gfortran -O2 Polzeros_src/pzeros.f90 ../src/fpml.f90 test_software.f90 -L AMVW_src -lAMVW -o test_software

# clean up
rm *.mod
rm *.o
rm *.a