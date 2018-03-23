#
# Makefile for FPML tests
# 
include ../../make.inc

all: conv init_est_acc init_est_time methods rand_poly spec_poly_list spec_poly_results unity clean
	
conv:
	$(LC) $(LFLAGS) conv.tex
	
init_est_acc:
	$(LC) $(LFLAGS) init_est_acc.tex
	
init_est_time:
	$(LC) $(LFLAGS) init_est_time.tex
	
methods:
	$(LC) $(LFLAGS) methods.tex
	
rand_poly:
	$(LC) $(LFLAGS) rand_poly.tex
	
spec_poly_list:
	$(LC) $(LFLAGS) spec_poly_list.tex
	
spec_poly_results:
	$(LC) $(LFLAGS) spec_poly_results.tex
	
unity:
	$(LC) $(LFLAGS) unity.tex
	
clean: copy
	rm *.log
	rm *.aux
	rm *.pdf
	rm *.png
	
copy:
	cp *.png ../figures/