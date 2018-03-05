#!/bin/bash

# compile fpml_driver
gfortran -O2 src/fpml.f90 src/fpml_driver.f90 -o fpml_driver

# clean up
rm *.mod