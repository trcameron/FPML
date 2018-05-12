#!/bin/bash

# export dyld library
export DYLD_LIBRARY_PATH=$PWD/AMVW/eiscor/lib:$DYLD_LIBRARY_PATH

# run tests
#echo "convergence test ..."
#./conv
#echo "initial estimates test ..."
#./init_est
#echo "methods test ... "
#./methods $1 $2 $3
echo "random polynomials test ..."
./rand_poly $1 $2 $3
echo "special polynomials test ..."
./spec_poly
echo "roots of unity test ..."
./unity $1 $2 $3