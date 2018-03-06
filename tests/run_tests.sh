#!/bin/bash

# run tests
echo "initial estimates test"
./init_est
echo "spec_poly test"
./spec_poly
echo "conv test"
./conv
echo "methods test, startDegree = 80, endDegree = 10240, maxit = 25"
./methods 80 10240 25
echo "rand_poly test, startDegree = 80, endDegree = 10240, maxit = 25"
./rand_poly 80 10240 25
echo "unity test, startDegree = 80, endDegree = 10240, maxit = 25"
./unity 80 10240 25