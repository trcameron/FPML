#!/bin/bash

# compile other_src
gfortran -c -O2 other_src/*.f90 -lblas
gfortran -c -O2 other_src/turnovers/*.f90 -lblas
rm poly_zeroes.mod
rm pzeros.o

# create library
ar crv libfpml_test.a *.o
cp libfpml_test.a lib/

# compile rand_poly.f90
gfortran -O2 other_src/pzeros.f90 ../src/fpml.f90 rand_poly.f90 -L lib -lfpml_test -lblas -o rand_poly
# compile spec_poly.f90
gfortran -O2 other_src/pzeros.f90 ../src/fpml.f90 spec_poly.f90 -L lib -lfpml_test -lblas -o spec_poly
# compile unity.f90
gfortran -O2 other_src/pzeros.f90 ../src/fpml.f90 unity.f90 -L lib -lfpml_test -lblas -o unity
# compile module.f90
gfortran -O2 ../src/fpml.f90 module.f90 -o module
# compile methods.f90
gfortran -O2 methods.f90 -o methods

# clean up
rm *.mod
rm *.o
rm *.a