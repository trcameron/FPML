#############################################
#           FPML Make File TeX Directory    #
#############################################
include ../../make.inc

all: conv nat_poly rand_poly rand_unity spec_poly_list spec_poly_results trunc_exp unity clean
	
conv:
	$(LC) $(LFLAGS) conv.tex
	
nat_poly:
	$(LC) $(LFLAGS) nat_poly.tex
	
rand_poly:
	$(LC) $(LFLAGS) rand_poly.tex
	
rand_unity:
	$(LC) $(LFLAGS) rand_unity.tex
	
spec_poly_list:
	$(LC) $(LFLAGS) spec_poly_list.tex
	
spec_poly_results:
	$(LC) $(LFLAGS) spec_poly_results.tex
	
trunc_exp:
	$(LC) $(LFLAGS) trunc_exp.tex
	
unity:
	$(LC) $(LFLAGS) unity.tex
	
clean: copy
	@rm -r *.log
	@rm -r *.aux
	@rm -r *.pdf
	@rm -r *.png
	
copy:
	@cp -r *.png ../figures/