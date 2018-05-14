#!/bin/sh

nagdir=/home/thomas/NAG/fll6i26dfl
nagldir="${nagdir}/lib"
nagmoddir="${nagdir}/nag_interface_blocks"
fcompile="gfortran -g -Wall -I${nagmoddir}"
flink="${nagldir}/libnag_nag.a -lstdc++"

$fcompile ../../src/fpml.f90 rand_poly_nag.f90 $flink -o rand_poly_nag.exe

$fcompile ../../src/fpml.f90 spec_poly_nag.f90 $flink -o spec_poly_nag.exe

$fcompile ../../src/fpml.f90 unity_nag.f90 $flink -o unity_nag.exe
