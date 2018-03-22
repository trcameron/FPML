#
# Makefile for FPML 
#
include make.inc

install:
	$(FC) $(FFLAGS) -o fpml_driver src/fpml.f90 src/fpml_driver.f90
	$(MAKE) all -C tests
	
uninstall: clean
	rm -f fpml_driver
	$(MAKE) uninstall -C tests

run:
	$(MAKE) run -C tests
	
clean:
	rm -f *.mod
	$(MAKE) clean -C tests