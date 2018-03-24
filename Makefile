#############################################
#           FPML Make File Root directory   #
#############################################
include make.inc

install:
	$(FC) $(FFLAGS) -o fpml_driver src/fpml.f90 src/fpml_driver.f90
	@$(MAKE) install -C tests
	
uninstall: clean
	@$(MAKE) uninstall -C tests
	@rm -f fpml_driver

run:
	@$(MAKE) run -C tests

compile:
	@$(MAKE) compile -C tests
	
clean:
	@$(MAKE) clean -C tests
	@rm -f *.mod