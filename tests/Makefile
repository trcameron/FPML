#############################################
#           FPML Make File Tests Directory  #
#############################################
include ../make.inc ../.master.inc

MPOBJS := mpfun/mpfun90.o mpfun/mpmod90.o mpfun/mpmodm90.o mpfun/mpmodx90.o

install: convg eiscor mpfunsrc nat_poly rand_poly rand_unity rev_nat_poly spec_poly trunc_exp unity
	
convg: mpfunsrc
	$(FC) $(FFLAGS) -o convg ../src/fpml.f90 src/rootstocoeffs.f90 src/convg.f90 $(MPOBJS) -I mpfun
	
nat_poly: eiscor
	$(FC) $(FFLAGS) -o nat_poly ../src/fpml.f90 Polzeros/pzeros.f90 src/nat_poly.f90 -I AMVW/include AMVW/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)

rand_poly: eiscor
	$(FC) $(FFLAGS) -o rand_poly ../src/fpml.f90 Polzeros/pzeros.f90 src/rand_poly.f90 -I AMVW/include AMVW/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)

rand_unity: eiscor mpfunsrc
	$(FC) $(FFLAGS) -o rand_unity ../src/fpml.f90 Polzeros/pzeros.f90 src/rootstocoeffs.f90 src/rand_unity.f90 $(MPOBJS) -I mpfun -I AMVW/include AMVW/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)

rev_nat_poly: eiscor
	$(FC) $(FFLAGS) -o rev_nat_poly ../src/fpml.f90 Polzeros/pzeros.f90 src/rev_nat_poly.f90 -I AMVW/include AMVW/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)

spec_poly: eiscor mpfunsrc
	$(FC) $(FFLAGS) -o spec_poly ../src/fpml.f90 Polzeros/pzeros.f90 src/rootstocoeffs.f90 src/spec_poly.f90 $(MPOBJS) -I mpfun -I AMVW/include AMVW/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)
	
trunc_exp: eiscor
	$(FC) $(FFLAGS) -o trunc_exp ../src/fpml.f90 Polzeros/pzeros.f90 src/trunc_exp.f90 -I AMVW/include AMVW/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)
	
unity: eiscor
	$(FC) $(FFLAGS) -o unity ../src/fpml.f90 Polzeros/pzeros.f90 src/unity.f90 -I AMVW/include AMVW/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)

eiscor:
	@$(MAKE) install -C AMVW
    
mpfunsrc: 
	@$(MAKE) all -C mpfun
	
run:
	@./run_tests.sh

compile:
	@$(MAKE) all -C TeX
	
uninstall: clean
	@$(MAKE) uninstall -C AMVW
	@rm -f conv
	@rm -f nat_poly
	@rm -f rand_poly
	@rm -f rand_unity
	@rm -f rev_nat_poly
	@rm -f spec_poly
	@rm -f trunc_exp
	@rm -f unity

clean:
	@$(MAKE) clean -C AMVW
	@$(MAKE) clean -C mpfun
	@rm -f *.mod