#!/bin/bash

# run tests
echo "module test"
./module
echo "methods test, startDegree = 100, endDegree = 6400, maxit = 25"
./methods
echo "rand_poly test, startDegree = 100, endDegree = 6400, maxit = 25"
./rand_poly 100 6400 25
echo "unity test, startDegree = 100, endDegree = 6400, maxit = 25"
./unity 100 6400 25
echo "spec_poly test"
./spec_poly