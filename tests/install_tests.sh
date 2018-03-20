#!/bin/bash

gfortran -O3 -cpp test_eiscor.f90 -o test_eiscor -I /Users/thcameron/Documents/eiscor/include /Users/thcameron/eiscor/lib/libeiscor.dylib.0.2.0


# compile other_src (except for poly_zeros which is included as module)
#gfortran -c -O2 other_src/*.f90 -lblas
#gfortran -c -O2 other_src/turnovers/*.f90 -lblas
#rm poly_zeroes.mod
#rm pzeros.o

# create library
#ar crv libfpml_test.a *.o
#cp libfpml_test.a lib/

# compile rand_poly.f90
#gfortran -O2 other_src/pzeros.f90 ../src/fpml.f90 rand_poly.f90 -L lib -lfpml_test -lblas -o rand_poly
# compile spec_poly.f90
#gfortran -O2 other_src/pzeros.f90 ../src/fpml.f90 spec_poly.f90 -L lib -lfpml_test -lblas -o spec_poly
# compile unity.f90
#gfortran -O2 other_src/pzeros.f90 ../src/fpml.f90 unity.f90 -L lib -lfpml_test -lblas -o unity
# compile init_est.f90
#gfortran -O2 other_src/pzeros.f90 ../src/fpml.f90 init_est.f90 -o init_est
# compile methods.f90
#gfortran -O2 methods.f90 -o methods
# compile conv.f90
#gfortran -O2 ../src/fpml.f90 conv.f90 -o conv

# clean up
rm *.mod
rm *.o
rm *.a