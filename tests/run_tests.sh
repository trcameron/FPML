#!/bin/bash

# export dyld library
export DYLD_LIBRARY_PATH=$PWD/AMVW/eiscor/lib:$DYLD_LIBRARY_PATH

# run tests
echo "convergence test ..."
./convg
#echo "natural polynomials test ..."
#./nat_poly
#echo "random polynomials test ..."
#./rand_poly
echo "random unity test ..."
./rand_unity
#echo "rev natural polynomials test ..."
#./rev_nat_poly
echo "special polynomials test ..."
./spec_poly
echo "truncated exponential test ..."
./trunc_exp
#echo "roots of unity test ..."
#./unity