#!/bin/bash

# run tests
echo "natural polynomials test ..."
./nat_poly_nag.exe
echo "random polynomials test ..."
./rand_poly_nag.exe
echo "random unity test ..."
./rand_unity_nag.exe
echo "rev natural polynomials test ..."
./rev_nat_poly_nag.exe
echo "special polynomials test ..."
./spec_poly_nag.exe
echo "truncated exponential test ..."
./trunc_exp_nag.exe
echo "roots of unity test ..."
./unity_nag.exe
