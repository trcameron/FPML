#############################################
#           FPML Make File TESTS DIRECTORY  #
#############################################
include ../make.inc

install: conv init_est methods rand_poly spec_poly unity eiscor
	
conv:
	$(FC) $(FFLAGS) -o conv ../src/fpml.f90 src/conv.f90
	
init_est:
	$(FC) $(FFLAGS) -o init_est ../src/fpml.f90 Polzeros/pzeros.f90 src/init_est.f90
	
methods:
	$(FC) $(FFLAGS) -o methods src/methods.f90
	
rand_poly: eiscor
	$(FC) $(FFLAGS) -o rand_poly ../src/fpml.f90 Polzeros/pzeros.f90 src/rand_poly.f90 -I AMVW/include AMVW/eiscor/lib/libeiscor.dylib.0.2.0

spec_poly: eiscor
	$(FC) $(FFLAGS) -o spec_poly ../src/fpml.f90 Polzeros/pzeros.f90 src/spec_poly.f90 -I AMVW/include AMVW/eiscor/lib/libeiscor.dylib.0.2.0
	
unity: eiscor
	$(FC) $(FFLAGS) -o unity ../src/fpml.f90 Polzeros/pzeros.f90 src/unity.f90 -I AMVW/include AMVW/eiscor/lib/libeiscor.dylib.0.2.0

eiscor:
	@$(MAKE) install -C AMVW
	
run:
	@./run_tests.sh $(SD) $(ED) $(IT)

compile:
	@$(MAKE) all -C TeX
	
uninstall: clean
	@$(MAKE) uninstall -C AMVW
	@rm -f conv
	@rm -f init_est
	@rm -f methods
	@rm -f rand_poly
	@rm -f spec_poly
	@rm -f unity

clean:
	@$(MAKE) clean -C AMVW
	@rm -f *.mod