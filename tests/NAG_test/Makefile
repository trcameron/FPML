#############################################
#           NAG_TEST MAKEFILE               #
#############################################
NAGDIR := /home/thomas/NAG/fll6i26dfl
NAGLDIR := $(NAGDIR)/lib
NAGMODDIR := $(NAGDIR)/nag_interface_blocks
FCOMPILE := gfortran -O2 -I $(NAGMODDIR)
FLINK := $(NAGLDIR)/libnag_nag.a -lstdc++
MPOBJS := mpfun/mpfun90.o mpfun/mpmod90.o mpfun/mpmodm90.o mpfun/mpmodx90.o

install: nat_poly rand_poly rand_unity rev_nat_poly spec_poly trunc_exp unity

nat_poly:
	$(FCOMPILE) -o nat_poly_nag.exe ../../src/fpml.f90 nat_poly_nag.f90 $(FLINK)

rand_poly:
	$(FCOMPILE) -o rand_poly_nag.exe ../../src/fpml.f90 rand_poly_nag.f90 $(FLINK)

rand_unity: mpfunsrc
	$(FCOMPILE) -o rand_unity_nag.exe ../../src/fpml.f90 rootstocoeffs.f90 rand_unity_nag.f90 $(FLINK) $(MPOBJS) -I mpfun

rev_nat_poly:
	$(FCOMPILE) -o rev_nat_poly_nag.exe ../../src/fpml.f90 rev_nat_poly_nag.f90 $(FLINK)

spec_poly: mpfunsrc
	$(FCOMPILE) -o spec_poly_nag.exe ../../src/fpml.f90 rootstocoeffs.f90 spec_poly_nag.f90 $(FLINK) $(MPOBJS) -I mpfun

trunc_exp:
	$(FCOMPILE) -o trunc_exp_nag.exe ../../src/fpml.f90 trunc_exp_nag.f90 $(FLINK)

unity:
	$(FCOMPILE) -o unity_nag.exe ../../src/fpml.f90 unity_nag.f90 $(FLINK)

mpfunsrc: 
	@$(MAKE) all -C mpfun

run:
	@./run_tests.sh

uninstall: clean
	@rm -f nat_poly_nag.exe
	@rm -f rand_poly_nag.exe
	@rm -f rand_unity_nag.exe
	@rm -f rev_nat_poly_nag.exe
	@rm -f spec_poly_nag.exe
	@rm -f trunc_exp_nag.exe
	@rm -f unity_nag.exe

clean:
	@$(MAKE) clean -C mpfun
	@rm -f *.mod
